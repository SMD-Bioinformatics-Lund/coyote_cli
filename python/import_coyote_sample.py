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
import re
from cli import cli_parser
from datetime import datetime
from bson.objectid import ObjectId

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
    update = False
    command = args.command_selection
    if command == "load":
        args_dict = vars(args)
        # remove None values from dict, otherwise it gets loaded into sample collection
        args_dict = {key: value for key, value in args_dict.items() if value is not None}
        # make groups a list. Load sample into more than one page in coyote. Coyote expects a list
        # although it is hardly ever used, and might cause error downstream with 2 or more groups
        tmp_list = []
        tmp_list.append(args_dict["groups"])
        args_dict["groups"] = tmp_list
    elif command == "yaml":
        args_dict = validate_yaml(args.yaml_file)
        args_dict["update"] = args.update
        args_dict["increment"] = args.increment
    sample_dict = {}
    # check what's being loaded, DNA or RNA
    data_type = data_typer(args_dict)
    for key in args_dict:
        if key in ["load","command_selection","debug_logger","quiet","increment","update"]:
            continue
        sample_dict[key] = args_dict[key]
    logging.debug(f"Sample meta information {sample_dict}")
    # do a load, get the ID-hash from sample load. Add this as SAMPLE_ID to all other documents per case
    client = pymongo.MongoClient(config.mongo["uri"])
    db = client[config.mongo["dbname"]]
    samples_col = db["samples"]
    # update case, get sample_id for sample-collection #
    if args_dict["update"]:
        logging.debug(f"Hi! Updating sample via {command} argument")
        update = True
        args_dict,sample_id = update_case(args_dict,samples_col,data_type)
        meta_info_updater(sample_dict,sample_id,samples_col)
    # NEW CASE check db for case_id, dont cause id crashes #
    else:
        logging.debug(f"Hi! \nLoading sample via {command} argument")
        sample_dict["name"] = what_id(args_dict["name"],args_dict["increment"],samples_col)
        sample_dict["time_added"] = datetime.utcnow()
        sample_id = samples_col.insert_one(sample_dict)
        sample_id = sample_id.inserted_id
    # Load DNA variation #
    if "vcf_files" in args_dict:
        # since this is a mandatory value for non-update cases. Need to be handled differently
        if args_dict["vcf_files"] != "no_update":
            logging.debug(f"Loading DNA variation, starting with SNV variants..")
            load_snvs(args_dict["vcf_files"],sample_id,args_dict["groups"],update,db)
            exit
        # load optional data
        if "cnv" in args_dict:
            logging.debug(f"Reading copy number variation")
            load_cnv_json(args_dict["cnv"],sample_id,update,db)
        if "biomarkers" in args_dict:
            logging.debug(f"Reading other biomarkers")
            load_biomarkers(args_dict["biomarkers"],sample_id,update,db)
        if "transloc" in args_dict:
            logging.debug(f"Reading DNA translocations")
            load_transloc(args_dict["transloc"],sample_id,update,db)
        if "lowcov" in args_dict:
            logging.debug(f"Reading regions with lower than expected coverage")
            load_lowcov(args_dict["lowcov"],sample_id,args_dict["name"],update,db)
    # Load RNA variation
    elif "fusion_files" in args_dict:
        # since this is a mandatory value for non-update cases. Need to be handled differently
        if args_dict["fusion_files"] != "no_update":
            logging.debug(f"Loading RNA variation, starting with fusions..")
            load_fusions(args_dict["fusion_files"],sample_id,update,db)
        # load optional data
        if "expression_path" in args_dict:
            logging.debug(f"Reading gene expression levels")
            with open(args_dict['expression_path'],'r') as file:
                exp = json.load(file)
                if update:
                    result = samples_col.update_one( {"_id": ObjectId(str(sample_id)) }, {'$unset': {'expr':""}})
                    logging.debug(f"Removed {result.modified_count} expression levels for sample")
                logging.debug(f"Inserted gene expression levels")
                samples_col.update_one( {"_id": ObjectId(str(sample_id)) }, {'$set': {'expr': exp}})
        if "classification_path" in args_dict:
            with open(args_dict['classification_path'],'r') as file:
                class_data = json.load(file)
                if update:
                    result = samples_col.update_one( {"_id": ObjectId(str(sample_id)) }, {'$unset': {'classification':""}})
                    logging.debug(f"Removed {result.modified_count} classifications for sample")
                samples_col.update_one( {"_id": ObjectId(str(sample_id)) }, {'$set': {'classification': class_data}})
                logging.debug(f"Reading classifications based upon expression levels")
        if "qc" in args_dict:
            logging.debug(f"Reading qc data")
            with open(args_dict["qc"],'r') as file:
                val = json.load(file)
                if update:
                    result = samples_col.update_one({"_id": ObjectId(str(sample_id))}, {"$unset": {"QC": ""}}) 
                    logging.debug(f"Removed {result.modified_count} QC for sample")
                samples_col.update_one( {"_id": ObjectId(str(sample_id))}, {"$set": {"QC": [val]}})
                logging.debug(f"Reading QC from the rnaseq pipeline")

def load_snvs(infile,sample_id,group,update,db):
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
        var_dict["SAMPLE_ID"] = str(sample_id)
        # find floats in CSQ, make sure they are saved into coyote as such
        # also, remove gnomad_AF&gnomad_AF types
        var_dict = emulate_perl(var_dict)
        # apply filters according to assay
        for key in config_filters:
            if "gnomAD" in key:
                for transcript in var_dict["INFO"]["CSQ"]:
                    # sometimes gnomAD AF is assigned as 0.01&0 use max AF of these
                    max_af = transcript["gnomAD_AF"]
                    if max_af == "":
                        max_af = 0.0
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
        # edits for coyote
        var_dict["INFO"]["variant_callers"] = var_dict["INFO"]["variant_callers"].split("|")
        var_dict["FILTER"] = var_dict["FILTER"].split(";")
        del var_dict["FORMAT"]
        count = 0
        for sample in var_dict["GT"]:
            if "AF" not in sample and "VAF" not in sample or "DP" not in sample or "VD" not in sample or "GT" not in sample:
                exit("not a valid VCF, should be aggregated by AF(VAF), VD AD and GT")
            # first sample is tumor, add this information to db, also change VAF to AF as coyote expects
            if not count:
                var_dict["GT"][count]["type"] = "case"
                var_dict["GT"][count]["AF"] =  var_dict["GT"][count]["VAF"]
                del var_dict["GT"][count]["VAF"]
            else:
                var_dict["GT"][count]["type"] = "control"
                var_dict["GT"][count]["AF"] =  var_dict["GT"][count]["VAF"]
                del var_dict["GT"][count]["VAF"]
            var_dict["GT"][count]["sample"] = var_dict["GT"][count]["_sample_id"]
            del var_dict["GT"][count]["_sample_id"]
            count += 1
        filtered_data.append(var_dict)
    if update:
        delete_collection("variants_idref",sample_id)
    collection = db["variants_idref"]
    result = collection.insert_many(filtered_data)
    logging.debug(f"{len(result.inserted_ids)} SNV/Indel variant imported")
    
def validate_yaml(yaml_file):
    """
    read YAML file. Confirm mandatory field are present and return dict akin to args
    """
    with open(yaml_file, 'r') as file:
        yaml_dict = yaml.safe_load(file)
    if ("vcf_files" not in yaml_dict or "fusion_files" not in yaml_dict) and "groups" not in yaml_dict and "name" not in yaml_dict and "genome_build" not in yaml_dict:
        exit("YAML is missing mandatory fields: vcf, groups, name or build")

    return yaml_dict

def load_cnv_json(cnv_json,sample_id,update,db):
    """
    read all variant in cnv json, add case_id to each variant
    save as list to import many to mongodb
    """
    cnv_variants = []
    with open(cnv_json,'r') as file:
        cnv_dict = json.load(file)
    for var in cnv_dict:
        cnv_dict[var]["SAMPLE_ID"] = str(sample_id)
        cnv_variants.append(cnv_dict[var])
    
    if update:
        delete_collection("cnvs_wgs",sample_id)
    collection = db["cnvs_wgs"]
    result = collection.insert_many(cnv_variants)
    logging.debug(f"Inserted {len(cnv_variants)} copy number variants")

def load_biomarkers(biomarkers_json,sample_id,update,db):
    """
    read biomarkers json, add SAMPLE_ID via case_id
    """
    with open(biomarkers_json,'r') as file:
        biomarkers_dict = json.load(file) 
    biomarkers_dict["SAMPLE_ID"] = str(sample_id)
    if update:
        delete_collection("biomarkers",sample_id)
    collection = db["biomarkers"]
    result = collection.insert_one(biomarkers_dict)
    logging.debug(f"Inserted {len(biomarkers_dict)-2} other biomarkers")

def load_lowcov(lowcov_bed,sample_id,case_id,update,db):
    """
    read lowcov bedfile and load to db with case_id as SAMPLE_ID
    """
    lowcov_data = []
    with open(lowcov_bed, 'r') as f:
        lowcov_dict = csv.DictReader(f,delimiter="\t",fieldnames=['chr', 'start', 'end', 'avg_cov', 'amplicon'])
        for row in lowcov_dict:
            row["SAMPLE_ID"] = str(sample_id)
            row["sample"] = case_id
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["avg_cov"] = float(row["avg_cov"])
            lowcov_data.append(row)
    if update:
        delete_collection("coverage",sample_id)
    collection = db["coverage"]
    result = collection.insert_many(lowcov_data)
    logging.debug(f"Inserted {len(lowcov_data)} regions with lower than expected coverage")

def load_transloc(infile,sample_id,update,db):
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
            var_dict["SAMPLE_ID"] = str(sample_id)
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
    if update:
        delete_collection("transloc",sample_id)
    collection = db["transloc"]
    result = collection.insert_many(filtered_data)
    logging.debug(f"Inserted {len(filtered_data)} DNA fusion variants")

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

def load_fusions(infile,sample_id,update,db):
    """
    read fusions JSON from rnaseq-fus pipeline
    add SAMPLE_ID to collections
    """
    fusions_list = []
    with open(infile,'r') as file:
        fusions_dict = json.load(file) 
    for fusion in fusions_dict:
        fusion["SAMPLE_ID"] = str(sample_id)
        fusions_list.append(fusion)
    if update:
        delete_collection("fusions",sample_id)
    collection = db["fusions"]
    result = collection.insert_many(fusions_list)
    logging.debug(f"Inserted {len(fusions_list)} RNA fusion variants")

def what_id(case_id,increment,samples_col):
    """
    lookup case id. If exist and increment=True, return new case_id with suffix
    if exist and increment=False, exit with status case exist
    if it's unique, return name again
    """
    samples_found_exact = list(samples_col.find( { "name": case_id } ))
    # exact name does not exist since before
    if not samples_found_exact:
        logging.debug(f"ID will be: {case_id}")
        return case_id
    if not increment:
        exit("Sample already exist, use --increment True if you want to add this id with autogenerated suffix, otherwise change ID")
    logging.debug(f"This id seems to be present with suffixes since before, trying to find a suitable suffix")
    suffixes = []
    true_matches = 0
    new_name = case_id
    add_suffix = ""
    samples_found = list(samples_col.find( { "name": { '$regex' : case_id }} ))
    # only accept matches with suffixes(right match), and not prefixes(left matches)
    for name in samples_found:
        left_match,right_match,true_match = catch_left_right(case_id,name["name"])
        if right_match and not left_match and true_match:
            suffixes.append(right_match)
            true_matches +=1
    max_suffix = 1
    if true_matches:
        # if for some reason there are multiple exact matches. Could exist old data in db sanity check
        if len(suffixes) == 0:
            exit("This sample seems to have beed added several times, consider fixing the db entry as this will cause error in displaying case in coyote")
        for suffix in suffixes:
            # if this id has been iterated before, find the max iteration. Ignore matches to id that has not been iterated.
            # if id-2 id-3 id-5 exist, next id will be id-6 and not id-4. id-2 id-3 exists new id will be id-4
            name_match = re.match('-\d+',suffix)
            if name_match:
                suffix_inter = int(suffix.replace("-",""))
                if suffix_inter > max_suffix:
                    max_suffix = suffix_inter
    add_suffix = str(max_suffix+1)
    new_name = f"{case_id}-{add_suffix}"
    logging.debug(f"New ID will be {new_name}")
    return new_name

def catch_left_right(case_id,name):
    """
    finds characters left and right of string match
    return right and left match
    """
    pattern = fr'(.*)({re.escape(case_id)})(.*)'
    matches = re.match(pattern, name)
    if matches:
        left_match = matches.group(1)
        true_match = matches.group(2)
        right_match = matches.group(3)
    return left_match,right_match,true_match

def update_case(args_dict,samples_col,data_type):
    """
    verify case exists, get objectID from sample-col
    verify that no DNA data is added to RNA sample and vice verse
    since fusion_files and vcf_files are mandatory for RNA respectively DNA
    these two are given no_update if they are not part of the update, but since
    they are mandatory they need a value
    could also add only meta-data, some of this is data agnostic, some are not
    """
    sample_id = ""
    try:
        samples_found_exact = dict(samples_col.find_one( { "name": args_dict["name"] } )) #, 'build':args_dict['build'] would make it impossible to update build
    except:
        exit("Cannot find case in database, will not update anything. Bye!")
    if samples_found_exact:
        sample_id = samples_found_exact["_id"]
        # find if it is RNA or DNA case
        if "fusion_files" in samples_found_exact:
            if data_type == "DNA":
                exit("you are trying to add DNA data to a RNA sample. BAD PERSON!")
            if "fusion_files" in args_dict:
                pass
            else:
                args_dict["fusion_files"] = "no_update"
        elif "vcf_files" in samples_found_exact:
            if data_type == "RNA":
                exit("you are trying to add RNA data to a DNA sample. BAD PERSON!")
            elif "vcf_files" in args_dict:
                pass
            else:
                args_dict["vcf_files"] = "no_update"
        else:
            exit("update function could not determine if the case is DNA or RNA")       
    return(args_dict,sample_id)

def delete_collection(collection_name,sample_id):
    """
    Deletes all entries for a given collection and matching SAMPLE_ID
    mathcing samples-collection ObjectID
    This is to facilitate update function
    """
    client = pymongo.MongoClient(config.mongo["uri"])
    db = client[config.mongo["dbname"]]
    collection = db[collection_name]
    sample_col = db["samples"]
    sample = dict(sample_col.find_one( {"_id": ObjectId(str(sample_id))} ))
    results = list(collection.find({ "SAMPLE_ID": str(sample_id) }))
    if results:
        results = collection.delete_many( {"SAMPLE_ID": str(sample_id)})
        logging.debug(f"{results.deleted_count} entries deleted from {collection_name} for {sample['name']}({sample_id})")
    else:
        logging.debug(f"No entries found for {collection_name} for {sample['name']}({sample_id})")

def data_typer(args_dict):
    """
    verify that one data type is being added to sample
    return type if successfull otherwise exit
    """
    data_type = None
    for dtype in args_dict:
        if dtype in config.data_types:
            if data_type is None:
                data_type = config.data_types[dtype]
            else:
                if config.data_types[dtype] != data_type:
                    exit("data types are both from RNA and DNA. Check your input")
    return data_type

def meta_info_updater(meta_dict,sample_id,samples_col):
    """
    check for changes, replace if new
    Add new field that did not exist before
    """
    try:
        result = dict(samples_col.find_one( {"_id": ObjectId(str(sample_id))}))
    except:
        exit("Cannot find case in database, will not update anything. Bye! meta info")
    for arg in meta_dict:
        if arg in result:
            if meta_dict[arg] != result[arg]:
                if arg == "group":
                    exit("No support to update group as of yet")
                samples_col.update_one( { "_id": ObjectId(str(sample_id)) }, { '$set': {str(arg): meta_dict[arg] }})
                logging.debug(f"changing {arg} : from {result[arg]} to {meta_dict[arg]}")
        else:
            samples_col.update_one( { "_id": ObjectId(str(sample_id)) }, { '$set': {str(arg): meta_dict[arg] }})
            logging.debug(f"adding {arg} : {meta_dict[arg]} to sample")

def emulate_perl(var_dict):
    """
    Perl is superior when it comes to auto detecting data types, let us emulate that
    """
    for transcript in var_dict["INFO"]["CSQ"]:
        for key in transcript:
            ## first try an split on &, and then first item for data type
            if isinstance(transcript[key], str):
                data = transcript[key].split("&")
                if is_float(data[0]):
                    data = [float(x) for x in data]
                    max_float = max(data)
                    transcript[key] = float(max_float)      
    return var_dict

def is_float(s):
    try:
        float(s)
        if len(s.split('.')) > 1:
            return True
    except ValueError:
        return False

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
        
        