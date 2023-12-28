from pysam import VariantFile
from pprint import pprint
import cmdvcf
import sys
import pymongo
import json
import logging
import argparse
import config
import yaml
import csv
import gzip
from cli import cli_parser

## Import Coyote Sample 
#  Adds sample to coyote database. Should be able to load the following data types
#  SNVs: Via VCF format, CSQ strings. Paired samples has to have tumor -> normal order in header
#  CNVs: Via JSON defining start end chrom log2 ratios and overlapping genes
#  Sample information: Via JSON, identifiers such as ID clarity ID, both for case and control
#  Translocations: Via VCF, SNPeff strings
#  Check whether sample exist, try to iterate new suffix
#  Low Cov data: Via Bed file
#  Biomarkers: Via JSON
#  Fusions: Via JSON
##

def main(args) -> None:
    command = args.command_selection
    print(command)
    if command == "load":
        args_dict = vars(args)
        args_dict = {key: value for key, value in args_dict.items() if value is not None}
        tmp_list = []
        tmp_list.append(args_dict["groups"])
        args_dict["groups"] = tmp_list
    elif command == "yaml" :
        args_dict = validate_yaml(args.yaml_file)
    # do a load, get the ID-hash from sample load. Add this as SAMPLE_ID to all other documents per case
    sample_id = "PH"
    sample_dict = {}
    for key in args_dict:
        if key in ["load","command_selection","debug_logger","quiet","increment"]:
            continue
        sample_dict[key] = args_dict[key]
    # client = pymongo.MongoClient(config.mongo["uri"])
    # db = client[config.mongo["dbname"]]
    # samples_col = db["samples"]
    sample_id = samples_col.insert_one(sample_dict)
    # get SNVs to load
    filtered_snvs = load_snvs(args_dict["vcf"],sample_id,args_dict["groups"])
    # load optional data
    if "cnv" in args_dict:
        cnv_variants = load_cnv_json(args_dict["cnv"],sample_id)
    # if "fusions" in args_dict:
    #     load_cnv_json()
    if "biomarkers" in args_dict:
        biomarkers = load_biomarkers(args_dict["biomarkers"],sample_id)
    if "transloc" in args_dict:
        load_transloc(args_dict["transloc"],sample_id)
    if "lowcov" in args_dict:
        load_lowcov(args_dict["lowcov"],sample_id,args_dict["name"])

def load_snvs(infile,sample_id,group):
    """
    A function to load variants into variants_idref. Only load usable information from CSQ!
    Path to VCF, will load cmdvcf module using pysam. In config per group(assay) define filters
    and what CSQ-fields to load.
    """
    filtered_data = []
    vcf_object = VariantFile(infile)
    for var in vcf_object.fetch():       
        var_dict = cmdvcf.parse_variant(var,vcf_object.header)
        filters = var_dict["FILTER"]
        # fix for pindel variants, add TYPE
        if "SVTYPE" in var_dict["INFO"]:
            var_dict["INFO"]["TYPE"] = var_dict["INFO"]["SVTYPE"]
        # filter according to config
        if group[0] in config.assay:
            config_filters = config.assay[group[0]]["filters"]
        else:
            config_filters = config.assay["default"]["filters"]
        fail_filter = 0
        # for DB match to present variants for this case
        var_dict["SAMPLE_ID"] = sample_id
        for key in config_filters:
            if "gnomAD" in key:
                for transcript in var_dict["INFO"]["CSQ"]:
                    # sometimes gnomAD AF is assigned as 0.01&0 use max AF of these
                    gnomadaf = transcript["gnomAD_AF"].split("&")
                    max_af = max(gnomadaf)
                    if float(max_af) > config_filters[key]:
                        fail_filter = 1
            elif key in filters:
                fail_filter = 1
                break
        if fail_filter:
            break
        # delete unused CSQ-fields using config
        count = 0
        for transcript in var_dict["INFO"]["CSQ"]:
            for field in config.coyote_csq:
                try:
                    del var_dict["INFO"]["CSQ"][count][field]
                except:
                    continue
            count +=1
        count = 0
        for sample in var_dict["GT"]:
            if "AF" not in sample and "VAF" not in sample or "DP" not in sample or "VD" not in sample or "GT" not in sample:
                print("not a valid VCF, should be aggregated by AF(VAF), VD AD and GT")
            # first sample is tumor, add this information to db
            if not count:
                var_dict["GT"][count]["type"] = "case"
            else:
                var_dict["GT"][count]["type"] = "control"
            count += 1
        filtered_data.append[var_dict]
    return filtered_data
    
def validate_yaml(yaml_file):
    """
    read YAML file. Confirm mandatory field are present and return dict akin to args
    """
    with open(yaml_file, 'r') as file:
        yaml_dict = yaml.safe_load(file)
    if "vcf" not in yaml_dict or "groups" not in yaml_dict or "name" not in yaml_dict or "build" not in yaml_dict:
        exit("YAML is missing mandatory fields: vcf, groups, name or build")

    return yaml_dict

def load_cnv_json(cnv_json,sample_id):
    """
    read all variant in cnv json, add case_id to each variant
    save as list to import many to mongodb
    """
    cnv_variants = []
    with open(cnv_json,'r') as file:
        cnv_dict = json.load(file)
    for var in cnv_dict:
        cnv_dict[var]["SAMPLE_ID"] = sample_id
        cnv_variants.append(cnv_dict[var])
    return cnv_variants

def load_biomarkers(biomarkers_json,sample_id):
    """
    read biomarkers json, add SAMPLE_ID via case_id
    """
    with open(biomarkers_json,'r') as file:
        biomarkers_dict = json.load(file) 
    biomarkers_dict["SAMPLE_ID"] = sample_id
    return biomarkers_dict

def load_lowcov(lowcov_bed,sample_id,case_id):
    """
    read lowcov bedfile and load to db with case_id as SAMPLE_ID
    """
    lowcov_data = []
    with open(lowcov_bed, 'r') as f:
        lowcov_dict = csv.DictReader(f,delimiter="\t",fieldnames=['chr', 'start', 'end', 'avg_cov', 'amplicon'])
        for row in lowcov_dict:
            row["SAMPLE_ID"] = sample_id
            row["sample"] = case_id
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["avg_cov"] = float(row["avg_cov"])
            lowcov_data.append(row)
    return lowcov_data

def load_transloc(infile,sample_id):
    """
    read VCF containing DNA fusions from manta-like BNDs
    Need to be annotated by SNPeff
    Add SAMPLE_ID hash key to each variant and return list of variants
    """
    mane = read_mane(config.mane)
    filtered_data = []
    vcf_object = VariantFile(infile)
    for var in vcf_object.fetch():
        var_dict = cmdvcf.parse_variant(var,vcf_object.header)
        # ignore dups and dels (must be a better way legacy from bjhall)
        if "<" not in var_dict["ALT"]:
            var_dict["SAMPLE_ID"] = sample_id
            keep_variant = 0
            mane_select = {}
            all_new_ann = []
            add_mane = 0
            for ann in var_dict["INFO"]["ANN"]:
                ## count mane matches for both genes in pair
                n_mane = 0
                genes = ann["Gene_ID"].split('&')
                for gene in genes:
                    enst = mane.get(gene, {}).get('ensembl', 'NO_MANE_TRANSCRIPT')
                    if enst in ann["HGVS.p"]:
                        n_mane +=1
                new_ann = {}
                ## keep bidirectional and fusion annotations
                for key in ann:
                    if key == "Annotation":
                        for anno in ann["Annotation"]:
                            if anno == "gene_fusion" or anno == "bidirectional_gene_fusion":
                                keep_variant = 1
                    ## a lot of dot notation in SNPeff, remove from final import to DB
                    dotless_key = key.replace(".","")
                    new_ann[dotless_key] = ann[key]
                all_new_ann.append(new_ann)
                # if both genes in pair are mane save for MANE_ANN annotation
                if n_mane > 0 and n_mane == len(genes):
                    mane_select = new_ann
                    add_mane = 1 
        del var_dict["INFO"]["ANN"]
        var_dict["INFO"]["ANN"] = all_new_ann
        if add_mane:
            var_dict["INFO"]["MANE_ANN"] = mane_select
        if keep_variant:
            filtered_data.append(var_dict)
    return filtered_data

def read_mane(txt_gz):
    """
    read gzipped tab-delimited file, create a mane to refseq and mane to ensembl library
    """
    mane_dict = {}
    with gzip.open(txt_gz, 'rt') as f:
        mane_file = csv.DictReader(f,delimiter="\t")
        for line in mane_file:
            refseq = line["RefSeq_nuc"].split(".")[0]
            ensembl = line["Ensembl_nuc"].split(".")[0]
            tmp = {}
            tmp["refseq"] = refseq
            tmp["ensembl"] = ensembl
            ensembl_gene = line["Ensembl_Gene"].split(".")[0]
            mane_dict[ensembl_gene] = tmp
    return mane_dict


def setup_logging(debug : bool = False) -> None:

    format = '[%(asctime)s][%(levelname)s]: %(message)s'
    if debug:
        logging.basicConfig(level = logging.DEBUG, format = format)
    else:
        logging.basicConfig(level = logging.INFO, format = format)


if __name__ == '__main__':

    parser = cli_parser()
    args = parser.parse_args()
    if not args.command_selection:
        parser.print_help()
        sys.exit(2)

    setup_logging(debug = args.debug_logger)


    main(args)
        
        