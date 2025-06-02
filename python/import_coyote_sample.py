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
#  Low Cov data: Via Bed file / or via JSON cov.
#  Biomarkers: Via JSON
#  Fusions: Via JSON
##


def main(args) -> None:
    update = False
    command = args.command_selection
    if command == "load":
        args_dict = vars(args)
        # remove None values from dict, otherwise it gets loaded into sample collection
        args_dict = {
            key: value for key, value in args_dict.items() if value is not None
        }
        # make groups a list. Load sample into more than one page in coyote. Coyote expects a list
        # although it is hardly ever used, and might cause error downstream with 2 or more groups
        tmp_list = []
        tmp_list.append(args_dict["groups"])
        args_dict["groups"] = tmp_list
        control_id = args_dict.get("control_id", "null")
        args_dict["control_id"] = control_id if control_id != "null" else None
    elif command == "yaml":
        args_dict = validate_yaml(args.yaml_file)
        control_id = args_dict.get("control_id", "null")
        args_dict["update"] = args.update
        args_dict["increment"] = args.increment
        args_dict["control_id"] = control_id if control_id != "null" else None
    sample_dict = {}
    # check what's being loaded, DNA or RNA
    data_type = data_typer(args_dict)
    for key in args_dict:
        if key in [
            "load",
            "command_selection",
            "debug_logger",
            "quiet",
            "increment",
            "update",
        ]:
            continue
        sample_dict[key] = args_dict[key]
    logging.debug(f"Sample meta information {sample_dict}")
    # do a load, get the ID-hash from sample load. Add this as SAMPLE_ID to all other documents per case
    client = pymongo.MongoClient(config.mongo["uri"])
    db = client[config.mongo["dbname"]]
    samples_col = db["samples"]
    canonical_col = db["refseq_canonical"]
    # update case, get sample_id for sample-collection #
    if args_dict["update"]:
        logging.debug(f"Hi! Updating sample via {command} argument")
        update = True
        args_dict, sample_id = update_case(args_dict, samples_col, data_type)
        meta_info_updater(sample_dict, sample_id, samples_col)
    # NEW CASE check db for case_id, dont cause id crashes #
    else:
        logging.debug(f"Hi! \nLoading sample via {command} argument")
        sample_dict["name"] = what_id(
            args_dict["name"], args_dict["increment"], samples_col
        )
        sample_dict["time_added"] = datetime.utcnow()
        sample_id = samples_col.insert_one(sample_dict)
        sample_id = sample_id.inserted_id
    # Load DNA variation #
    if "vcf_files" in args_dict:
        # since this is a mandatory value for non-update cases. Need to be handled differently
        if args_dict["vcf_files"] != "no_update":
            logging.debug(f"Loading DNA variation, starting with SNV variants..")
            load_snvs(
                args_dict["vcf_files"], sample_id, args_dict["groups"], update, db
            )
            exit
        # load optional data
        if "cnv" in args_dict:
            logging.debug(f"Reading copy number variation")
            load_cnv_json(args_dict["cnv"], sample_id, update, db)
        if "biomarkers" in args_dict:
            logging.debug(f"Reading other biomarkers")
            load_biomarkers(args_dict["biomarkers"], sample_id, update, db)
        if "transloc" in args_dict:
            logging.debug(f"Reading DNA translocations")
            load_transloc(args_dict["transloc"], sample_id, update, db)
        if "lowcov" in args_dict and "cov" in args_dict:
            logging.debug(f"Both lowcov and cov is being loaded for sample, pick one")
            exit()
        if "lowcov" in args_dict:
            logging.debug(f"Reading regions with lower than expected coverage")
            load_lowcov(args_dict["lowcov"], sample_id, args_dict["name"], update, db)
        if "cov" in args_dict:
            logging.debug(f"Reading coverage data from JSON")
            load_cov(args_dict["cov"], sample_id, args_dict["name"], update, db)
    # Load RNA variation
    elif "fusion_files" in args_dict:
        # since this is a mandatory value for non-update cases. Need to be handled differently
        if args_dict["fusion_files"] != "no_update":
            logging.debug(f"Loading RNA variation, starting with fusions..")
            load_fusions(args_dict["fusion_files"], sample_id, update, db)
        # load optional data
        if "expression_path" in args_dict:
            logging.debug(f"Reading gene expression levels")
            with open(args_dict["expression_path"], "r") as file:
                exp = json.load(file)
                if update:
                    result = samples_col.update_one(
                        {"_id": ObjectId(str(sample_id))}, {"$unset": {"expr": ""}}
                    )
                    logging.debug(
                        f"Removed {result.modified_count} expression levels for sample"
                    )
                logging.debug(f"Inserted gene expression levels")
                samples_col.update_one(
                    {"_id": ObjectId(str(sample_id))}, {"$set": {"expr": exp}}
                )
        if "classification_path" in args_dict:
            with open(args_dict["classification_path"], "r") as file:
                class_data = json.load(file)
                if update:
                    result = samples_col.update_one(
                        {"_id": ObjectId(str(sample_id))},
                        {"$unset": {"classification": ""}},
                    )
                    logging.debug(
                        f"Removed {result.modified_count} classifications for sample"
                    )
                samples_col.update_one(
                    {"_id": ObjectId(str(sample_id))},
                    {"$set": {"classification": class_data}},
                )
                logging.debug(f"Reading classifications based upon expression levels")
        if "qc" in args_dict:
            logging.debug(f"Reading qc data")
            with open(args_dict["qc"], "r") as file:
                val = json.load(file)
                if update:
                    result = samples_col.update_one(
                        {"_id": ObjectId(str(sample_id))}, {"$unset": {"QC": ""}}
                    )
                    logging.debug(f"Removed {result.modified_count} QC for sample")
                samples_col.update_one(
                    {"_id": ObjectId(str(sample_id))}, {"$set": {"QC": [val]}}
                )
                logging.debug(f"Reading QC from the rnaseq pipeline")


def refseq_noversion(acc):
    a = acc.split(".")
    return a[0]


def select_csq(csq_arr, canonical):

    db_canonical = -1
    vep_canonical = -1
    first_protcoding = -1

    impact_order = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

    for impact in impact_order:
        for csq_idx, csq in enumerate(csq_arr):
            if csq["IMPACT"] == impact:

                if csq["SYMBOL"] in canonical and canonical[
                    csq["SYMBOL"]
                ] == refseq_noversion(csq["Feature"]):
                    db_canonical = csq_idx
                    return (csq_arr[db_canonical], "db")
                if csq["CANONICAL"] == "YES" and vep_canonical == -1:
                    vep_canonical = csq_idx
                if (
                    first_protcoding == -1
                    and csq["BIOTYPE"] == "protein_coding"
                    and first_protcoding == -1
                ):
                    first_protcoding = csq_idx

    if vep_canonical >= 0:
        return (csq_arr[vep_canonical], "vep")
    elif first_protcoding >= 0:
        return (csq_arr[first_protcoding], "random")

    return (csq_arr[0], "random")


def selcted_transcript_removal(csq_arr, selected_transcript):
    """
    remove all seected transcript from the other transcripts
    """
    for csq_index, csq in enumerate(csq_arr):
        if csq["Feature"] == selected_transcript:
            del csq_arr[csq_index]
            break
    return csq_arr


## CANONICAL TRANSCRIPTS
def get_canonical():
    canonical_dict = {}
    client = pymongo.MongoClient(config.mongo["uri"])
    db = client[config.mongo["dbname"]]
    canonical_col = db["refseq_canonical"]
    canonical = canonical_col.find({})
    for c in canonical:
        canonical_dict[c["gene"]] = c["canonical"]

    return canonical_dict


canonical_dict = get_canonical()


def load_snvs(infile, sample_id, group, update, db):
    """
    A function to load variants into variants_idref. Only load usable information from CSQ!
    Path to VCF, will load cmdvcf module using pysam. In config per group(assay) define filters
    and what CSQ-fields to load.
    """
    filtered_data = []
    vcf_object = VariantFile(infile)
    for var in vcf_object.fetch():
        var_dict = cmdvcf.parse_variant(var, vcf_object.header)
        var_csq = var_dict["INFO"]["CSQ"]
        # removes the variant all together if the csq only has experimental transcripts but retains if the any of those experimental transcripts are in the list of genes
        if var_csq:
            all_features = [c.get("Feature") for c in var_csq]
            all_X_genes = [
                c.get("SYMBOL") for c in var_csq if c.get("Feature", "").startswith("X")
            ]

        if all([f.startswith("X") for f in all_features]) and not any(
            [
                g in ["HNF1A", "MZT2A", "SNX9", "KLHDC4", "LMTK3", "PTPA"]
                for g in list(set(all_X_genes))
            ]
        ):
            continue

        # fix for pindel variants, add TYPE
        if "SVTYPE" in var_dict["INFO"]:
            var_dict["INFO"]["TYPE"] = var_dict["INFO"]["SVTYPE"]
        # for DB match to present variants for this case
        var_dict["SAMPLE_ID"] = str(sample_id)
        # find floats in CSQ, make sure they are saved into coyote as such
        var_dict = emulate_perl(var_dict)
        # edits for coyote
        # pick AF field to present in collection.
        var_dict.update(pick_af_fields(var_dict))
        var_dict["variant_class"] = var_csq[0].get("VARIANT_CLASS")
        # summerize variant for easier indexing and searching between annot-collection and variants_idref
        (
            slim_csq,
            cosmic_list,
            dbsnp,
            pubmed_list,
            transcripts_list,
            cdna_list,
            prot_list,
            genes_list,
            hotspots_dict,
        ) = parse_transcripts(var_csq)

        # select csq based on old logic
        selected_csq, selected_csq_source = select_csq(
            slim_csq, canonical_dict
        )  # How slow is this?

        var_dict["INFO"]["CSQ"] = selcted_transcript_removal(
            slim_csq, selected_csq["Feature"]
        )
        var_dict["INFO"]["selected_CSQ"] = selected_csq
        var_dict["INFO"]["selected_CSQ_criteria"] = selected_csq_source
        var_dict["selected_csq_feature"] = selected_csq["Feature"]
        var_dict["HGVSp"] = prot_list
        var_dict["HGVSc"] = cdna_list
        var_dict["genes"] = genes_list
        var_dict["transcripts"] = transcripts_list
        var_dict["INFO"]["CSQ"] = slim_csq
        var_dict["cosmic_ids"] = cosmic_list
        var_dict["dbsnp_id"] = dbsnp
        var_dict["pubmed_ids"] = pubmed_list
        var_dict["hotspots"] = [hotspots_dict]
        var_dict["simple_id"] = (
            f"{var_dict['CHROM']}_{var_dict['POS']}_{var_dict['REF']}_{var_dict['ALT']}"
        )
        var_dict["INFO"]["variant_callers"] = var_dict["INFO"]["variant_callers"].split(
            "|"
        )
        var_dict["FILTER"] = var_dict["FILTER"].split(";")
        del var_dict["FORMAT"]
        count = 0
        for sample in var_dict["GT"]:
            if (
                "AF" not in sample
                and "VAF" not in sample
                or "DP" not in sample
                or "VD" not in sample
                or "GT" not in sample
            ):
                exit("not a valid VCF, should be aggregated by AF(VAF), VD AD and GT")
            # first sample is tumor, add this information to db, also change VAF to AF as coyote expects
            if not count:
                var_dict["GT"][count]["type"] = "case"
                var_dict["GT"][count]["AF"] = var_dict["GT"][count]["VAF"]
                del var_dict["GT"][count]["VAF"]
            else:
                var_dict["GT"][count]["type"] = "control"
                var_dict["GT"][count]["AF"] = var_dict["GT"][count]["VAF"]
                del var_dict["GT"][count]["VAF"]
            var_dict["GT"][count]["sample"] = var_dict["GT"][count]["_sample_id"]
            del var_dict["GT"][count]["_sample_id"]
            count += 1
        filtered_data.append(var_dict)
    if update:
        delete_collection("variants_idref", sample_id)
    collection = db["variants"]
    result = collection.insert_many(filtered_data)
    logging.debug(f"{len(result.inserted_ids)} SNV/Indel variant imported")


def validate_yaml(yaml_file):
    """
    read YAML file. Confirm mandatory field are present and return dict akin to args
    """
    with open(yaml_file, "r") as file:
        yaml_dict = yaml.safe_load(file)
    if (
        ("vcf_files" not in yaml_dict or "fusion_files" not in yaml_dict)
        and "groups" not in yaml_dict
        and "name" not in yaml_dict
        and "genome_build" not in yaml_dict
    ):
        exit("YAML is missing mandatory fields: vcf, groups, name or build")

    return yaml_dict


def load_cnv_json(cnv_json, sample_id, update, db):
    """
    read all variant in cnv json, add case_id to each variant
    save as list to import many to mongodb
    """
    cnv_variants = []
    with open(cnv_json, "r") as file:
        cnv_dict = json.load(file)
    for var in cnv_dict:
        cnv_dict[var]["SAMPLE_ID"] = str(sample_id)
        cnv_variants.append(cnv_dict[var])

    if update:
        delete_collection("cnvs_wgs", sample_id)
    collection = db["cnvs_wgs"]
    if len(cnv_variants) > 0:
        result = collection.insert_many(cnv_variants)
        logging.debug(f"Inserted {len(cnv_variants)} copy number variants")
    else:
        logging.debug(f"No CNVs to insert")


def load_biomarkers(biomarkers_json, sample_id, update, db):
    """
    read biomarkers json, add SAMPLE_ID via case_id
    """
    with open(biomarkers_json, "r") as file:
        biomarkers_dict = json.load(file)
    biomarkers_dict["SAMPLE_ID"] = str(sample_id)
    if update:
        delete_collection("biomarkers", sample_id)
    collection = db["biomarkers"]
    result = collection.insert_one(biomarkers_dict)
    logging.debug(f"Inserted {len(biomarkers_dict)-2} other biomarkers")


def load_lowcov(lowcov_bed, sample_id, case_id, update, db):
    """
    read lowcov bedfile and load to db with case_id as SAMPLE_ID
    """
    lowcov_data = []
    with open(lowcov_bed, "r") as f:
        lowcov_dict = csv.DictReader(
            f, delimiter="\t", fieldnames=["chr", "start", "end", "avg_cov", "amplicon"]
        )
        for row in lowcov_dict:
            row["SAMPLE_ID"] = str(sample_id)
            row["sample"] = case_id
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["avg_cov"] = float(row["avg_cov"])
            lowcov_data.append(row)
    if update:
        delete_collection("coverage", sample_id)
    collection = db["coverage"]
    result = collection.insert_many(lowcov_data)
    logging.debug(
        f"Inserted {len(lowcov_data)} regions with lower than expected coverage"
    )


def load_cov(cov_json, sample_id, case_id, update, db):
    """
    read coverage JSON-file and load to db with case_id as SAMPLE_ID
    """
    with open(cov_json, "r") as file:
        cov_dict = json.load(file)
    cov_dict["SAMPLE_ID"] = str(sample_id)
    cov_dict["sample"] = str(case_id)
    if update:
        delete_collection("panel_cov", sample_id)
    collection = db["panel_cov"]
    result = collection.insert_one(cov_dict)
    logging.debug(f"Inserted coverage data")


def load_transloc(infile, sample_id, update, db):
    """
    read VCF containing DNA fusions from manta-like BNDs
    Need to be annotated by SNPeff
    Add SAMPLE_ID hash key to each variant and return list of variants
    """
    mane = read_mane(config.mane)
    filtered_data = []
    vcf_object = VariantFile(infile)
    for var in vcf_object.fetch():
        var_dict = cmdvcf.parse_variant(var, vcf_object.header)
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
                genes = ann["Gene_ID"].split("&")
                for gene in genes:
                    enst = mane.get(gene, {}).get("ensembl", "NO_MANE_TRANSCRIPT")
                    if enst in ann["HGVS.p"]:
                        n_mane += 1
                new_ann = {}
                ## keep bidirectional and fusion annotations
                for key in ann:
                    if key == "Annotation":
                        for anno in ann["Annotation"]:
                            if (
                                anno == "gene_fusion"
                                or anno == "bidirectional_gene_fusion"
                            ):
                                keep_variant = 1
                    ## a lot of dot notation in SNPeff, remove from final import to DB
                    dotless_key = key.replace(".", "")
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
        delete_collection("transloc", sample_id)
    collection = db["transloc"]
    if len(filtered_data) > 0:
        result = collection.insert_many(filtered_data)
        logging.debug(f"Inserted {len(filtered_data)} DNA fusion variants")
    else:
        logging.debug(f"No DNA fusion variants to insert")


def read_mane(txt_gz):
    """
    read gzipped tab-delimited file, create a mane to refseq and mane to ensembl library
    """
    mane_dict = {}
    with gzip.open(txt_gz, "rt") as f:
        mane_file = csv.DictReader(f, delimiter="\t")
        for line in mane_file:
            refseq = line["RefSeq_nuc"].split(".")[0]
            ensembl = line["Ensembl_nuc"].split(".")[0]
            tmp = {}
            tmp["refseq"] = refseq
            tmp["ensembl"] = ensembl
            ensembl_gene = line["Ensembl_Gene"].split(".")[0]
            mane_dict[ensembl_gene] = tmp
    return mane_dict


def load_fusions(infile, sample_id, update, db):
    """
    read fusions JSON from rnaseq-fus pipeline
    add SAMPLE_ID to collections
    """
    fusions_list = []
    with open(infile, "r") as file:
        fusions_dict = json.load(file)
    for fusion in fusions_dict:
        fusion["SAMPLE_ID"] = str(sample_id)
        fusions_list.append(fusion)
    if update:
        delete_collection("fusions", sample_id)
    collection = db["fusions"]
    result = collection.insert_many(fusions_list)
    logging.debug(f"Inserted {len(fusions_list)} RNA fusion variants")


def what_id(case_id, increment, samples_col):
    """
    lookup case id. If exist and increment=True, return new case_id with suffix
    if exist and increment=False, exit with status case exist
    if it's unique, return name again
    """
    samples_found_exact = list(samples_col.find({"name": case_id}))
    # exact name does not exist since before
    if not samples_found_exact:
        logging.debug(f"ID will be: {case_id}")
        return case_id
    if not increment:
        exit(
            "Sample already exist, use --increment True if you want to add this id with autogenerated suffix, otherwise change ID"
        )
    logging.debug(
        f"This id seems to be present with suffixes since before, trying to find a suitable suffix"
    )
    suffixes = []
    true_matches = 0
    new_name = case_id
    add_suffix = ""
    samples_found = list(samples_col.find({"name": {"$regex": case_id}}))
    # only accept matches with suffixes(right match), and not prefixes(left matches)
    for name in samples_found:
        left_match, right_match, true_match = catch_left_right(case_id, name["name"])
        if right_match and not left_match and true_match:
            suffixes.append(right_match)
            true_matches += 1
    max_suffix = 1
    if true_matches:
        # if for some reason there are multiple exact matches. Could exist old data in db sanity check
        if len(suffixes) == 0:
            exit(
                "This sample seems to have beed added several times, consider fixing the db entry as this will cause error in displaying case in coyote"
            )
        for suffix in suffixes:
            # if this id has been iterated before, find the max iteration. Ignore matches to id that has not been iterated.
            # if id-2 id-3 id-5 exist, next id will be id-6 and not id-4. id-2 id-3 exists new id will be id-4
            name_match = re.match("-\d+", suffix)
            if name_match:
                suffix_inter = int(suffix.replace("-", ""))
                if suffix_inter > max_suffix:
                    max_suffix = suffix_inter
    add_suffix = str(max_suffix + 1)
    new_name = f"{case_id}-{add_suffix}"
    logging.debug(f"New ID will be {new_name}")
    return new_name


def catch_left_right(case_id, name):
    """
    finds characters left and right of string match
    return right and left match
    """
    pattern = rf"(.*)({re.escape(case_id)})(.*)"
    matches = re.match(pattern, name)
    if matches:
        left_match = matches.group(1)
        true_match = matches.group(2)
        right_match = matches.group(3)
    return left_match, right_match, true_match


def update_case(args_dict, samples_col, data_type):
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
        samples_found_exact = dict(
            samples_col.find_one({"name": args_dict["name"]})
        )  # , 'build':args_dict['build'] would make it impossible to update build
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
    return (args_dict, sample_id)


def delete_collection(collection_name, sample_id):
    """
    Deletes all entries for a given collection and matching SAMPLE_ID
    mathcing samples-collection ObjectID
    This is to facilitate update function
    """
    client = pymongo.MongoClient(config.mongo["uri"])
    db = client[config.mongo["dbname"]]
    collection = db[collection_name]
    sample_col = db["samples"]
    sample = dict(sample_col.find_one({"_id": ObjectId(str(sample_id))}))
    results = list(collection.find({"SAMPLE_ID": str(sample_id)}))
    if results:
        results = collection.delete_many({"SAMPLE_ID": str(sample_id)})
        logging.debug(
            f"{results.deleted_count} entries deleted from {collection_name} for {sample['name']}({sample_id})"
        )
    else:
        logging.debug(
            f"No entries found for {collection_name} for {sample['name']}({sample_id})"
        )


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


def meta_info_updater(meta_dict, sample_id, samples_col):
    """
    check for changes, replace if new
    Add new field that did not exist before
    """
    try:
        result = dict(samples_col.find_one({"_id": ObjectId(str(sample_id))}))
    except:
        exit("Cannot find case in database, will not update anything. Bye! meta info")
    for arg in meta_dict:
        if arg in result:
            if meta_dict[arg] != result[arg]:
                if arg == "group":
                    exit("No support to update group as of yet")
                samples_col.update_one(
                    {"_id": ObjectId(str(sample_id))},
                    {"$set": {str(arg): meta_dict[arg]}},
                )
                logging.debug(
                    f"changing {arg} : from {result[arg]} to {meta_dict[arg]}"
                )
        else:
            samples_col.update_one(
                {"_id": ObjectId(str(sample_id))}, {"$set": {str(arg): meta_dict[arg]}}
            )
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


def pick_af_fields(var):
    """
    return an AF field to top-level for variant

    Will have gnomAD_AF then gnomAD genomes, then exac, then thousand genomes
    and if possible return max af for gnomad
    """
    af_dict = {
        "gnomad_frequency": "",
        "gnomad_max": "",
        "exac_frequency": "",
        "thousandG_frequency": "",
    }
    allele = var["ALT"]
    exac = parse_allele_freq(var["INFO"]["CSQ"][0].get("ExAC_MAF"), allele)
    thousand_g = parse_allele_freq(var["INFO"]["CSQ"][0].get("GMAF"), allele)
    gnomad = var["INFO"]["CSQ"][0].get("gnomAD_AF", 0)
    gnomad_genome = var["INFO"]["CSQ"][0].get("gnomADg_AF", 0)
    gnomad_max = var["INFO"]["CSQ"][0].get("MAX_AF", 0)

    if gnomad:
        gnomad = max_gnomad(gnomad)
        af_dict["gnomad_frequency"] = gnomad
        if gnomad_max:
            af_dict["gnomad_max"] = gnomad_max
    elif gnomad_genome:
        gnomad = max_gnomad(gnomad)
        af_dict["gnomad_frequency"] = gnomad_genome
        if gnomad_max:
            af_dict["gnomad_max"] = gnomad_max
    else:
        af_dict["gnomad_frequency"] = ""
        af_dict["gnomad_max"] = ""
    if exac:
        af_dict["exac_frequency"] = exac
    if thousand_g:
        af_dict["thousandG_frequency"] = thousand_g
    return af_dict


def max_gnomad(gnomad):
    """
    check if gnoamd is multivalued, split and max
    """
    try:
        gnomad_list = gnomad.split("&")
        if gnomad_list:
            return float(max(gnomad_list))
    except:
        return gnomad


def parse_allele_freq(freq_str, allele):

    if freq_str:
        all_alleles = freq_str.split("&")
        for allele_frq in all_alleles:
            a = allele_frq.split(":")
            if a[0] == allele:
                return float(a[1])

    return 0


def split_on_comma(data):

    data_parts = data.split(":")
    if len(data_parts) > 1:
        return data_parts[1]
    else:
        return data


def split_on_ambersand(found_dict, string):
    """
    collect pubmed ids
    """
    try:
        string_list = string.split("&")
        for c in string_list:
            found_dict[c] = 1
        return found_dict
    except:
        found_dict[str(string)] = 1
        return found_dict


def collect_dbsnp(dbsnp_dict, dbsnp):
    """
    collect all rsids
    """
    dbsnp_list = dbsnp.split("&")
    for snp in dbsnp_list:
        if snp.startswith("rs"):
            dbsnp_dict[snp] = 1
    return dbsnp_dict


def collect_hotspots(hotspot_dict):
    cleaned_hotspot_dict = {}
    for hotspot, ids in hotspot_dict.items():
        formatted_ids = list(set(filter(None, ids)))
        if formatted_ids:
            cleaned_hotspot_dict[hotspot] = formatted_ids
    return cleaned_hotspot_dict


def parse_transcripts(csq: dict):
    """
    reduce redundancy from CSQ-strings
    re-organize data into top-level
    """
    transcripts = []
    pubmed_dict = {}
    cosmic_dict = {}
    dbsnp_dict = {}
    transcript_ids = {}
    hgvsc_ids = {}
    hgvsp_ids = {}
    gene_symbols = {}
    hotspots = {}

    for transcript in csq:
        slim_transcript = {}

        slim_transcript["Feature"] = transcript.get("Feature")
        transcript_id = str(transcript.get("Feature")).split(".")[0]
        transcript_ids[transcript_id] = 1
        slim_transcript["HGNC_ID"] = transcript.get("HGNC_ID")
        gene_symbol = transcript.get("SYMBOL")
        slim_transcript["SYMBOL"] = gene_symbol
        gene_symbols[gene_symbol] = 1
        slim_transcript["PolyPhen"] = transcript.get("PolyPhen")
        slim_transcript["SIFT"] = transcript.get("SIFT")
        slim_transcript["Consequence"] = transcript.get("Consequence")
        slim_transcript["ENSP"] = transcript.get("ENSP")
        slim_transcript["BIOTYPE"] = transcript.get("BIOTYPE")
        slim_transcript["INTRON"] = transcript.get("INTRON")
        slim_transcript["EXON"] = transcript.get("EXON")
        slim_transcript["CANONICAL"] = transcript.get("CANONICAL")
        slim_transcript["MANE"] = transcript.get("MANE_SELECT")
        slim_transcript["STRAND"] = transcript.get("STRAND")
        slim_transcript["IMPACT"] = transcript.get("IMPACT")
        slim_transcript["CADD_PHRED"] = transcript.get("CADD_PHRED")
        slim_transcript["CLIN_SIG"] = transcript.get("CLIN_SIG")
        slim_transcript["VARIANT_CLASS"] = transcript.get("VARIANT_CLASS")

        protein_change = transcript.get("HGVSp")
        if protein_change:
            protein_change = split_on_comma(protein_change)
        slim_transcript["HGVSp"] = protein_change
        hgvsp_ids[protein_change] = 1
        cdna_change = transcript.get("HGVSc")
        if cdna_change:
            cdna_change = split_on_comma(cdna_change)
        slim_transcript["HGVSc"] = cdna_change
        hgvsc_ids[cdna_change] = 1
        cosmic = transcript.get("COSMIC")
        if cosmic:
            cosmic_dict = split_on_ambersand(cosmic_dict, cosmic)

        dbsnp = transcript.get("Existing_variation")
        if dbsnp:
            dbsnp_dict = collect_dbsnp(dbsnp_dict, dbsnp)

        pubmed = transcript.get("PUBMED")
        if pubmed:
            pubmed_dict = split_on_ambersand(pubmed_dict, pubmed)

        for transcript_key in transcript.keys():
            for hsp in ["d", "gi", "lu", "cns", "mm", "co"]:
                if f"{hsp}hotspot_OID" in transcript_key:
                    hotspot = transcript.get(transcript_key)
                    if hotspot:
                        hotspots.setdefault(hsp, []).append(hotspot)

        transcripts.append(slim_transcript)

    cosmic_list = list(cosmic_dict.keys())
    pubmed_list = list(pubmed_dict.keys())
    dbsnp = list(dbsnp_dict.keys())
    hotspot_oids = collect_hotspots(hotspots)

    ## summerized
    transcript_list = list(transcript_ids.keys())
    transcript_list_filtered = [item for item in transcript_list if item]

    hgvsc_list = list(hgvsc_ids.keys())
    hgvsc_list_filtered = [item for item in hgvsc_list if item]

    hgvsp_list = list(hgvsp_ids.keys())
    hgvsp_list_filtered = [item for item in hgvsp_list if item]

    genes_list = list(gene_symbols.keys())
    genes_list_filtered = [item for item in genes_list if item]
    if len(dbsnp) > 0:
        dbsnp = dbsnp[0]
    else:
        dbsnp = ""

    return (
        transcripts,
        cosmic_list,
        dbsnp,
        pubmed_list,
        transcript_list_filtered,
        hgvsc_list_filtered,
        hgvsp_list_filtered,
        genes_list_filtered,
        hotspot_oids,
    )


def is_float(s):
    try:
        float(s)
        if len(s.split(".")) > 1:
            return True
    except ValueError:
        return False


def setup_logging(debug: bool = False) -> None:

    format = "[%(asctime)s][%(levelname)s]: %(message)s"
    if debug:
        logging.basicConfig(level=logging.DEBUG, format=format)
    else:
        logging.basicConfig(level=logging.INFO, format=format)


if __name__ == "__main__":

    parser = cli_parser()
    args = parser.parse_args()
    if not args.command_selection:
        parser.print_help()
        sys.exit(2)

    setup_logging(debug=args.debug_logger)

    main(args)
