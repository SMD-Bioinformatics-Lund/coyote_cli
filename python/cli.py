import argparse
import sys


MISSING_DOCUMENTATION = (
    "[HELP MISSING] Help text for this command/option is missing. Help out by making a "
    "pull request at: github.com/clinical-Genomics-Lund/bnf-infrastructure/"
)


def cli_parser():

    parser = argparse.ArgumentParser(
        description="Use/add/modify the data in the coyote database!",
        epilog="Show help for each subcommand with import_coyote_sample.py {subcommand} -h",
    )

    subparsers = parser.add_subparsers(dest="command_selection")
    base_subparser = argparse.ArgumentParser(add_help=False)

    """
    GLOBAL COMMANDS
    """
    base_subparser.add_argument(
        "--debug",
        action="store_true",
        dest="debug_logger",
        default=False,
        help="Print debug information to STDOUT.",
    )

    base_subparser.add_argument(
        "-q",
        "--quiet",
        default=False,
        action="store_true",
        help="When loading via CRON jobs hide output",
    )

    """
    MAIN COMMANDS
    """

    usage_msg = "For usage examples, see: http://mtlucmds1.lund.skane.se/wiki/doku.php?id=coyote"

    load = subparsers.add_parser(
        "load",
        parents=[base_subparser],
        description="Load case into coyote, flag for each datatype",
        help="Load case into coyote, flag for each datatype. EITHER RNA or DNA data types, never both!",
        epilog=usage_msg,
    )
    yaml = subparsers.add_parser(
        "yaml",
        parents=[base_subparser],
        description="Load case into coyote, all input defined in yaml-format",
        help="Load case into coyote, all input defined in yaml-format",
        epilog=usage_msg,
    )
    """
    YAML Custom args
    """
    yaml.add_argument(
        "-y",
        "--yaml",
        dest="yaml_file",
        required=True,
        help="Path to YAML file with all sample meta info and relevant VCFs JSONs and BED files to load",
    )
    yaml.add_argument(
        "--increment",
        dest="increment",
        default=True,
        help="If case ID already exists, increment this case by number of matches + 1",
    )
    yaml.add_argument(
        "-u",
        "--update",
        dest="update",
        default=False,
        action="store_true",
        help="Update existing case with new information or add new variation type",
    )
    """
    LOAD Custom args
    """
    sample_meta = load.add_argument_group("Sample Meta Information")
    dna_data = load.add_argument_group("Paths to DNA data")
    rna_data = load.add_argument_group("Paths to RNA data")
    ## Sample Meta Information ##
    sample_meta.add_argument(
        "--id",
        dest="name",
        required=True,
        help="Case ID, name of case in coyote. Related flags --increment (required)",
    )
    sample_meta.add_argument(
        "--group",
        dest="groups",
        required=True,
        help="Group name of case, to which group of coyote does sample belong. Example: myeloid_GMSv1.0 (required)",
    )
    sample_meta.add_argument(
        "--clarity-sample-id",
        dest="clarity-sample-id",
        help="Clarity Sample ID, need to make clarity db matches",
    )
    sample_meta.add_argument(
        "--clarity-pool-id",
        dest="clarity-pool-id",
        help="Of which pool was this sample part of",
    )
    sample_meta.add_argument(
        "--gens",
        dest="gens_id",
        help="ID of this case in GENS cmdscout2.lund.skane.se/gens",
    )
    sample_meta.add_argument(
        "--subpanel",
        dest="subpanel",
        help="Is this sample subcategorized into a smaller panel? Example for solid_GMSv3.0: ovarian|solid|breast|cns|lungthyroid",
    )
    sample_meta.add_argument(
        "--purity",
        dest="purity",
        type=float,
        help="Tumour purity of sample",
    )
    sample_meta.add_argument(
        "--build",
        dest="genome_build",
        required=True,
        help="To which genome build is the case aligned to (required)",
    )
    sample_meta.add_argument(
        "--increment",
        dest="increment",
        default=True,
        help="If case ID already exists, increment this case by number of matches + 1",
    )
    sample_meta.add_argument(
        "--tumor_sample",
        dest="tumor_sample",
        help="tumor_sample name",
    )
    sample_meta.add_argument(
        "--normal_sample",
        dest="normal_sample",
        help="normal_sample name, if available",
    )
    sample_meta.add_argument(
        "--profile",
        dest="profile",
        help="environment profile to use, e.g. 'production', 'testing', 'development', 'validation'",
    )
    sample_meta.add_argument(
        "--sample_no",
        dest="sample_no",
        type=int,
        help="Sample number, user to indicated if it is a paired sample or not",
    )
    sample_meta.add_argument(
        "--assay",
        dest="assay",
        help="Assay name, e.g. 'solid_GMSv3'",
    )
    ## Paths to DNA data ##
    dna_data.add_argument(
        "--vcf",
        dest="vcf_files",
        help="Path to SNV VCF, VEP annotated. Preferred SomaticPanelPipeline annotation for maximum usage",
    )
    dna_data.add_argument(
        "--cnv",
        dest="cnv",
        help="Path to CNV JSON, required output from SomaticPanelPipeline",
    )
    dna_data.add_argument(
        "--transloc",
        dest="transloc",
        help="Path to SNPeff annotated DNA translocations",
    )
    dna_data.add_argument(
        "--lowcov",
        dest="lowcov",
        help="Path to annotated BED file with average read depths",
    )
    dna_data.add_argument(
        "--cov",
        dest="lowcov",
        help="Path coverage JSON file with average read depths",
    )
    dna_data.add_argument(
        "--biomarkers",
        dest="biomarkers",
        help="Path to biomarkers JSON",
    )
    dna_data.add_argument(
        "--cnvprofile",
        dest="cnvprofile",
        help="Path to cnvprofile image",
    )
    ## Path to RNA data ##
    rna_data.add_argument(
        "--fusions",
        dest="fusion_files",
        help="Path to JSON containing gene fusions generated by nextflow_rnaseqfusion",
    )
    rna_data.add_argument(
        "--expression",
        dest="expression_path",
        help="Path to JSON containing gene normalized gene expression levels",
    )
    rna_data.add_argument(
        "--classification",
        dest="classification_path",
        help="Path to JSON containing gene expression classifications",
    )
    rna_data.add_argument(
        "--qc",
        dest="qc",
        help="Path to JSON qc data metrices for the pipeline",
    )
    load.add_argument(
        "-u",
        "--update",
        dest="update",
        default=False,
        action="store_true",
        help="Update existing case with new information or add new variation type",
    )
    # GetOptions( \%opt, 'vcf=s', 'id=s', 'clarity-sample-id=s', 'clarity-pool-id=s', 'bam=s', 'group=s', 'cnv=s', 'transloc=s', 'qc=s', 'cnvprofile=s', 'build=s', 'gens=s', 'biomarkers=s', 'subpanel=s', 'purity=s', 'lowcov=s' );

    return parser
