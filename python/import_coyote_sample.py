"""
===========================================================================
Class-based Atomic Sample Ingestion
===========================================================================
"""

from __future__ import annotations

import os
import re
import csv
import gzip
import json
import logging
import pymongo
import sys
import yaml
from pysam import VariantFile
from argparse import ArgumentParser, Namespace
from typing import Dict, Any, Generator, Optional, Tuple, Protocol, Literal
from dataclasses import dataclass
from contextlib import contextmanager
from datetime import datetime
from bson.objectid import ObjectId

# External project deps expected:
# - config (mongo URIs, db names, data_types mapping, mane path)
# - cmdvcf  (variant parser)
# - cli.cli_parser (arg parser)
import config
import cmdvcf
from cli import cli_parser


# --------------------------
# Logging
# --------------------------
def setup_logging(debug: bool = False) -> None:
    fmt = "[%(asctime)s][%(levelname)s]: %(message)s"
    logging.basicConfig(level=(logging.DEBUG if debug else logging.INFO), format=fmt)


# --------------------------
# Repositories (DB access)
# --------------------------
@dataclass
class Repos:
    """
    A class to manage interactions with a MongoDB database.
    Attributes:
        client (pymongo.MongoClient): The MongoDB client instance.
        db (pymongo.database.Database): The MongoDB database instance.
    """

    client: pymongo.MongoClient
    db: pymongo.database.Database

    @classmethod
    def from_args(cls, args_dict: Dict[str, Any]) -> "Repos":
        """
        Creates an instance of the Repos class using the provided arguments.

        Args:
            args_dict (Dict[str, Any]): A dictionary containing arguments.
            The key "dev" determines whether to use the development database
            or the production database.

        Returns:
            Repos: An instance of the Repos class initialized with a MongoDB client
            and the appropriate database.

        Notes:
            - The MongoDB client is created using the URI specified in the
            configuration (config.mongo["uri"]).
            - The database name is determined based on the "dev" key in args_dict:
            - If "dev" is True, the development database name (config.mongo["dbname_dev"]) is used.
            - Otherwise, the production database name (config.mongo["dbname"]) is used.
        """
        client = pymongo.MongoClient(config.mongo["uri"])
        dbname = (
            config.mongo["dbname_dev"]
            if args_dict.get("dev")
            else config.mongo["dbname"]
        )
        logging.info(f"Using database: {dbname}")
        return cls(client=client, db=client[dbname])

    @property
    def samples(self) -> pymongo.collection.Collection:
        return self.db["samples"]

    def col(self, name: str) -> pymongo.collection.Collection:
        return self.db[name]


# --------------------------
# Utilities
# --------------------------
CASE_CONTROL_KEYS: list[str] = [
    "case_id",
    "control_id",
    "clarity_control_id",
    "clarity_case_id",
    "clarity_case_pool_id",
    "clarity_control_pool_id",
    "case_ffpe",
    "control_ffpe",
    "case_sequencing_run",
    "control_sequencing_run",
    "case_reads",
    "control_reads",
    "case_purity",
    "control_purity",
]


def normalize_case_control(
    args: Dict[str, Any],
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Extracts and normalizes case and control sub-dictionaries from a flat dictionary of arguments.

    This function processes a dictionary of arguments (`args`) and separates it into two
    dictionaries: one for "case" keys and another for "control" keys. Keys starting with
    "case_" are added to the `case` dictionary (with the "case_" prefix removed), and keys
    starting with "control_" are added to the `ctrl` dictionary (with the "control_" prefix removed).

    Additionally, if any of the keys in `CASE_CONTROL_KEYS` are present in `args` and their
    value is either `None` or the string `"null"`, their value is normalized to `None`.

    Args:
        args (Dict[str, Any]): A flat dictionary of arguments containing keys for both
            "case" and "control" data.

    Returns:
        Tuple[Dict[str, Any], Dict[str, Any]]: A tuple containing two dictionaries:
            - `case`: A dictionary with normalized "case" data.
            - `ctrl`: A dictionary with normalized "control" data.
    """
    for k in CASE_CONTROL_KEYS:
        if k in args:
            if args[k] is None or args[k] == "null":
                args[k] = None
    case: dict = {}
    ctrl: dict = {}
    for k in CASE_CONTROL_KEYS:
        if "case" in k:
            case[k.replace("case_", "")] = args.get(k)
        elif "control" in k:
            ctrl[k.replace("control_", "")] = args.get(k)
    return case, ctrl


def build_sample_meta_dict(args_dict: Dict[str, Any]) -> Dict[str, Any]:
    """
    Constructs a metadata dictionary for a sample based on the provided arguments.
    This function processes the input dictionary `args_dict` to create a structured
    dictionary containing metadata for a sample. It separates case and control data
    using the `normalize_case_control` function and excludes certain keys from the
    resulting dictionary based on predefined rules.
    Args:
        args_dict (Dict[str, Any]): A dictionary containing input arguments and metadata.
    Returns:
        Dict[str, Any]: A dictionary containing the structured metadata for the sample.
                        Includes keys for "case" and optionally "control" if a control ID
                        is provided in the input.
    """

    sample_dict: Dict[str, Any] = {}
    case_dict, control_dict = normalize_case_control(args_dict)

    for key, value in args_dict.items():
        if key in [
            "load",
            "command_selection",
            "debug_logger",
            "quiet",
            "increment",
            "update",
            "dev",
        ]:
            continue
        if key in CASE_CONTROL_KEYS and key not in ["case_id", "control_id"]:
            continue
        sample_dict[key] = value

    sample_dict["case"] = case_dict
    if args_dict.get("control_id"):
        sample_dict["control"] = control_dict
    return sample_dict


def data_typer(args_dict: Dict[str, Any]) -> Optional[str]:
    """
    Determines the data type based on the keys in the provided dictionary.
    This function iterates through the keys in the input dictionary and checks if they
    match any of the predefined data types in the `config.data_types` mapping. If multiple
    data types are found and they conflict (e.g., both RNA and DNA), the program exits with
    an error message. If a single data type is identified, it is returned.
    Args:
        args_dict (Dict[str, Any]): A dictionary where the keys are checked against
            predefined data types in `config.data_types`.
    Returns:
        Optional[str]: The identified data type if found, otherwise `None`.
    Raises:
        SystemExit: If conflicting data types (e.g., RNA and DNA) are found in the input.
    """
    dtype = None
    for k in args_dict:
        if k in config.data_types:
            if dtype is None:
                dtype: str = config.data_types[k]
            elif config.data_types[k] != dtype:
                exit(
                    "Data types conflict: both RNA and DNA detected. Check your input."
                )
    return dtype


def validate_yaml(yaml_file: str) -> Dict[str, Any]:
    """
    Validate and parse a YAML file.

    This function reads a YAML file, checks for the presence of mandatory fields,
    and returns its contents as a dictionary. If any of the mandatory fields
    ("vcf_files", "fusion_files", "groups", "name", or "genome_build") are missing,
    the program will terminate with an error message.

    Args:
        yaml_file (str): The path to the YAML file to be validated.

    Returns:
        Dict[str, Any]: A dictionary containing the parsed contents of the YAML file.

    Raises:
        SystemExit: If the YAML file is missing any mandatory fields.
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


def catch_left_right(case_id: str, name: str) -> Tuple[str, str, str]:
    """
    Splits the input string `name` into three parts: the portion before the `case_id`,
    the `case_id` itself, and the portion after the `case_id`.

    Args:
        case_id (str): The identifier to search for within the `name` string.
        name (str): The input string to be split.

    Returns:
        Tuple[str, str, str]: A tuple containing three strings:
            - The portion of `name` before `case_id` (left_match).
            - The portion of `name` after `case_id` (right_match).
            - The `case_id` itself (true_match). If no match is found, all three strings will be empty.
    """
    left_match = right_match = true_match = ""
    pattern: str = rf"(.*)({re.escape(case_id)})(.*)"
    match: re.Match[str] | None = re.match(pattern, name)
    if match:
        left_match, true_match, right_match = (
            match.group(1),
            match.group(2),
            match.group(3),
        )
    return left_match, right_match, true_match


def what_id(case_id: str, increment: bool, samples_col) -> str:
    """
    Generate or validate a unique case ID based on the provided parameters.

    This function checks if the given `case_id` already exists in the `samples_col` collection.
    If it exists and `increment` is set to True, it generates a new unique ID by appending a numeric suffix.
    If it exists and `increment` is False, the function exits with an error message.
    If the `case_id` is unique, it returns the same `case_id`.

    Args:
        case_id (str): The case ID to validate or increment.
        increment (bool): Whether to generate a new ID with a suffix if the `case_id` already exists.
        samples_col: A database collection object used to query existing case IDs.

    Returns:
        str: A unique case ID, either the original or a new one with a numeric suffix.

    Raises:
        SystemExit: If `increment` is False and the `case_id` already exists, or if there are multiple exact matches
                    in the database that prevent generating a unique suffix.

    Notes:
        - The function assumes that `samples_col` supports MongoDB-like queries.
        - The suffix is determined by finding the highest numeric suffix among existing IDs and incrementing it.
        - If no suffixes are found, the function defaults to appending `-1` to the `case_id`.
    """
    existing_exact = list(samples_col.find({"name": case_id}))
    if not existing_exact:
        logging.info(f"Loading sample with the case ID: {case_id}")
        return case_id
    if not increment:
        exit("Sample already exist, use --increment True to auto-suffix, or change ID")

    logging.info("Existing suffixes present, choosing next suffix")
    suffixes: list = []
    true_matches: int = 0
    for doc in samples_col.find({"name": {"$regex": case_id}}):
        left, right, true = catch_left_right(case_id, doc["name"])
        if right and not left and true:
            suffixes.append(right)
            true_matches += 1

    max_suffix = 1
    if true_matches:
        if len(suffixes) == 0:
            exit("Multiple exact matches; fix DB to avoid display errors")
        for s in suffixes:
            m: re.Match[str] | None = re.match("-\d+", s)
            if m:
                n = int(s.replace("-", ""))
                if n > max_suffix:
                    max_suffix: int = n

    new_name = f"{case_id}-{max_suffix + 1}"
    logging.info(f"Generated new unique case ID: {new_name}")
    return new_name


def _exists(p: Optional[str]) -> bool:
    """
    Checks if the given path exists in the file system.

    Args:
        p (Optional[str]): The file or directory path to check. Can be None.

    Returns:
        bool: True if the path is not None and exists, False otherwise.
    """
    return bool(p) and os.path.exists(p)


def _require_exists(label: str, path: Optional[str]) -> None:
    """
    Ensures that the specified path exists; raises FileNotFoundError if not.

    Args:
        label (str): A descriptive label for the path, used in the error message.
        path (Optional[str]): The file or directory path to check. Can be None.

    Raises:
        FileNotFoundError: If the path is None or does not exist.
    """
    if not _exists(path):
        raise FileNotFoundError(f"{label} missing or not readable: {path}")


def _is_no_update(val: Optional[str]) -> bool:
    """
    Checks if the given value is a string that, when stripped of whitespace
    and converted to lowercase, matches the string "no_update".

    Args:
        val (Optional[str]): The value to check. It can be a string or None.

    Returns:
        bool: True if the value is a string equal to "no_update"
            (case-insensitive and ignoring surrounding whitespace),
            otherwise False.
    """
    return isinstance(val, str) and val.strip().lower() == "no_update"


def refseq_noversion(acc: str) -> str:
    """
    Removes the version suffix from a RefSeq accession number.

    A RefSeq accession number typically includes a version suffix
    separated by a dot (e.g., "NM_001256789.1"). This function
    extracts and returns only the base accession number without
    the version.

    Args:
        acc (str): The RefSeq accession number, which may include
                    a version suffix.

    Returns:
        str: The base accession number without the version suffix.

    Example:
        >>> refseq_noversion("NM_001256789.1")
        'NM_001256789'
    """
    return acc.split(".")[0]


def is_float(s: str) -> bool:
    """
    Determine if the given string represents a floating-point number.

    Args:
        s (str): The string to check.

    Returns:
        bool: True if the string represents a floating-point number (contains a decimal point and can be converted to a float), False otherwise.
    """
    try:
        float(s)
        return len(s.split(".")) > 1
    except ValueError:
        logging.debug(f"ValueError converting to float: {s}")
        return False
    except TypeError:
        logging.debug(f"TypeError converting to float: {s}")
        return False
    except Exception as e:
        logging.debug(f"Unexpected error converting to float: {s}, {e}")
        return False


def emulate_perl(var_dict) -> Any:
    """
    Processes and modifies the input dictionary by iterating over the "CSQ" transcripts
    in the "INFO" field. For each key in a transcript, if the value is a string, it splits
    the string on the '&' character and attempts to convert the first item to a float. If
    successful, it replaces the value with the maximum float value derived from the split
    data.

    Args:
        var_dict (dict): A dictionary containing an "INFO" field with a "CSQ" key, which
                        is expected to be a list of dictionaries representing transcripts.

    Returns:
        dict: The modified input dictionary with updated transcript values where applicable.

    Notes:
        - Assumes the presence of a helper function `is_float` that checks if a string can
        be converted to a float.
        - If the value for a key in a transcript is not a string, it is left unmodified.
    """
    for transcript in var_dict["INFO"]["CSQ"]:
        for key in list(transcript.keys()):
            ## first try an split on &, and then first item for data type
            if isinstance(transcript[key], str):
                data: Any = transcript[key].split("&")
                if is_float(data[0]):
                    transcript[key] = float(max(float(x) for x in data))
    return var_dict


def split_on_comma(data: Optional[str]) -> Optional[str]:
    """
    Splits the input string on the first occurrence of a colon (":")
    and returns the part after the colon. If the input string does not
    contain a colon or is empty/None, the original input is returned.

    Args:
        data (Optional[str]): The input string to be split. Can be None.

    Returns:
        Optional[str]: The part of the string after the first colon,
        or the original input if no colon is found or the input is None.
    """
    if not data:
        return data
    parts: list[str] = data.split(":")
    return parts[1] if len(parts) > 1 else data


def split_on_ambersand(found_dict: Dict[str, int], string) -> Dict[str, int]:
    """
    Splits the input string by the '&' character and updates the provided dictionary
    with each substring as a key and the value set to 1. If an exception occurs,
    the entire string is added to the dictionary as a key with the value set to 1.

    Args:
        found_dict (Dict[str, int]): A dictionary to store the substrings as keys
                                    with their values set to 1.
        string (str): The input string to be split by the '&' character.

    Returns:
        Dict[str, int]: The updated dictionary containing the substrings as keys
                        with their values set to 1.
    """
    try:
        string_list: Any = string.split("&")
        for c in string_list:
            found_dict[c] = 1
        return found_dict
    except Exception as e:
        logging.debug(f"Error processing string '{string}': {e}")
        found_dict[str(string)] = 1
        return found_dict


def collect_dbsnp(dbsnp_dict: Dict[str, int], dbsnp: str) -> Dict[str, int]:
    """
    Collects dbSNP identifiers from the input string and updates the provided dictionary.

    Args:
        dbsnp_dict (Dict[str, int]): A dictionary to store dbSNP identifiers.
        dbsnp (str): The input string containing dbSNP identifiers separated by '&'.

    Returns:
        Dict[str, int]: The updated dictionary containing dbSNP identifiers.
    """
    for snp in dbsnp.split("&"):
        if snp.startswith("rs"):
            dbsnp_dict[snp] = 1
    return dbsnp_dict


def collect_hotspots(hotspot_dict: Dict[str, list]) -> Dict[str, list]:
    """
    Splits the input hotspot dictionary and collects unique IDs for each hotspot.

    Args:
        hotspot_dict (Dict[str, list]): A dictionary containing hotspots and their associated IDs.

    Returns:
        Dict[str, list]: A dictionary with cleaned hotspot entries.
    """
    cleaned_hotspot_dict: dict = {}
    for hotspot, ids in hotspot_dict.items():
        formatted_ids: list = list(set(filter(None, ids)))
        if formatted_ids:
            cleaned_hotspot_dict[hotspot] = formatted_ids
    return cleaned_hotspot_dict


def parse_transcripts(csq: list) -> Any:
    """
    Parses and processes a list of transcript data, reducing redundancy and
    reorganizing the data into a structured format.

    Args:
        csq (list): A list of dictionaries containing transcript data. Each dictionary
                    represents a transcript with various attributes.

    Returns:
        tuple: A tuple containing the following elements:
            - transcripts (list): A list of processed transcript dictionaries with
            selected attributes.
            - cosmic_list (list): A list of unique COSMIC identifiers.
            - dbsnp_first (str): The first dbSNP identifier, if available.
            - pubmed_list (list): A list of unique PUBMED identifiers.
            - transcript_list (list): A list of unique transcript IDs.
            - hgvsc_list (list): A list of unique HGVSc identifiers.
            - hgvsp_list (list): A list of unique HGVSp identifiers.
            - gene_list (list): A list of unique gene symbols.
            - hotspot_oids (dict): A dictionary of hotspot OIDs categorized by type.

    Notes:
        - The function extracts and processes specific attributes from each transcript,
        such as "Feature", "HGNC_ID", "SYMBOL", "HGVSp", "HGVSc", and others.
        - It also collects and deduplicates identifiers from fields like "COSMIC",
        "PUBMED", and "Existing_variation".
        - Hotspot-related data is grouped by specific prefixes and returned as a dictionary.
    """
    transcripts: list = []
    pubmed_dict: dict = {}
    cosmic_dict: dict = {}
    dbsnp_dict: dict = {}
    transcript_ids: dict = {}
    hgvsc_ids: dict = {}
    hgvsp_ids: dict = {}
    gene_symbols: dict = {}
    hotspots: dict = {}

    for transcript in csq:
        slim: dict = {}
        feature: Any = transcript.get("Feature")
        slim["Feature"] = feature
        tid: str = str(feature).split(".")[0] if feature else ""
        if tid:
            transcript_ids[tid] = 1
        slim["HGNC_ID"] = transcript.get("HGNC_ID")
        sym: str = transcript.get("SYMBOL")
        slim["SYMBOL"] = sym
        if sym:
            gene_symbols[sym] = 1

        for k in (
            "PolyPhen",
            "SIFT",
            "Consequence",
            "ENSP",
            "BIOTYPE",
            "INTRON",
            "EXON",
            "CANONICAL",
            "MANE_SELECT",
            "STRAND",
            "IMPACT",
            "CADD_PHRED",
            "CLIN_SIG",
            "VARIANT_CLASS",
        ):
            v: Any = transcript.get(k)
            slim["MANE" if k == "MANE_SELECT" else k] = v

        protein: str | None = split_on_comma(transcript.get("HGVSp"))
        slim["HGVSp"] = protein
        if protein:
            hgvsp_ids[protein] = 1
        cdna: str | None = split_on_comma(transcript.get("HGVSc"))
        slim["HGVSc"] = cdna
        if cdna:
            hgvsc_ids[cdna] = 1

        cosmic: Any = transcript.get("COSMIC")
        if cosmic:
            cosmic_dict: Dict[str, int] = split_on_ambersand(cosmic_dict, cosmic)
        ev: Any = transcript.get("Existing_variation")
        if ev:
            dbsnp_dict: Dict[str, int] = collect_dbsnp(dbsnp_dict, ev)
        pm: Any = transcript.get("PUBMED")
        if pm:
            pubmed_dict: Dict[str, int] = split_on_ambersand(pubmed_dict, pm)

        for trk in list(transcript.keys()):
            for h in ["d", "gi", "lu", "cns", "mm", "co"]:
                if f"{h}hotspot_OID" in trk:
                    hv = transcript.get(trk)
                    if hv:
                        hotspots.setdefault(h, []).append(hv)

        transcripts.append(slim)

    cosmic_list = list(cosmic_dict.keys())
    pubmed_list = list(pubmed_dict.keys())
    dbsnp_list = list(dbsnp_dict.keys())
    dbsnp_first: str = dbsnp_list[0] if dbsnp_list else ""
    hotspot_oids: Dict[str, list] = collect_hotspots(hotspots)

    transcript_list: list = [t for t in transcript_ids.keys() if t]
    hgvsc_list: list = [t for t in hgvsc_ids.keys() if t]
    hgvsp_list: list = [t for t in hgvsp_ids.keys() if t]
    gene_list: list = [t for t in gene_symbols.keys() if t]

    return (
        transcripts,
        cosmic_list,
        dbsnp_first,
        pubmed_list,
        transcript_list,
        hgvsc_list,
        hgvsp_list,
        gene_list,
        hotspot_oids,
    )


def selected_transcript_removal(csq_arr, selected_transcript) -> Any:
    """
    Remove the selected transcript from the list of consequence annotations.

    Args:
        csq_arr (list): A list of consequence annotations (dictionaries
                            containing transcript information).
        selected_transcript (str): The transcript ID to be removed from the list.

    Returns:
        list: The updated list of consequence annotations with the selected
                transcript removed.
    """

    for idx, csq in enumerate(csq_arr):
        if csq.get("Feature") == selected_transcript:
            del csq_arr[idx]
            break
    return csq_arr


def select_csq(csq_arr, canonical: Dict[str, str]) -> Tuple[Dict[str, Any], str]:
    """
    Selects the most relevant consequence from a list of consequences based on impact
    and canonical status.

    The function prioritizes consequences in the following order:
    1. A consequence where the gene symbol matches the canonical dictionary and the
    feature matches the reference sequence without version.
    2. A consequence marked as "CANONICAL" by VEP (Variant Effect Predictor).
    3. The first consequence with a "protein_coding" biotype.
    4. The first consequence in the list if none of the above conditions are met.

    Args:
        csq_arr (list[Dict[str, Any]]): A list of dictionaries, where each dictionary
            represents a consequence with keys such as "IMPACT", "SYMBOL", "CANONICAL",
            "BIOTYPE", and "Feature".
        canonical (Dict[str, str]): A dictionary mapping gene symbols to their canonical
            reference sequences (without version).

    Returns:
        Tuple[Dict[str, Any], str]: A tuple containing the selected consequence dictionary
        and a string indicating the selection method:
            - "db": Selected based on the canonical dictionary.
            - "vep": Selected based on VEP's canonical annotation.
            - "random": Selected as the first protein-coding consequence or the first
            consequence in the list.
    """
    db_canonical = -1
    vep_canonical = -1
    first_protcoding = -1

    impact_order: list[str] = ["HIGH", "MODERATE", "LOW", "MODIFIER"]

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


def parse_allele_freq(freq_str, allele) -> float | Literal[0]:
    """
    Parses the allele frequency string for a specific allele.

    Args:
        freq_str (str): A string representing allele frequencies in the format
                        "A:0.01&C:0.02&T:0.03", where each allele and its frequency
                        are separated by '&' and each pair is separated by ':'.
        allele (str): The allele to search for in the frequency string.

    Returns:
        float: The frequency of the specified allele if found.
        Literal[0]: Returns 0 if the specified allele is not found or if the input
                    frequency string is empty.
    """
    if freq_str:
        for af in freq_str.split("&"):
            a: Any = af.split(":")
            if a[0] == allele:
                return float(a[1])
    return 0


def max_gnomad(gnomad) -> float | Any | None:
    """
    Processes a string containing values separated by '&', extracts the values,
    and returns the maximum value as a float. If an error occurs during processing,
    the original input is returned.

    Args:
        gnomad (str): A string containing values separated by '&'.

    Returns:
        float: The maximum value among the extracted parts, converted to a float.
        Any: The original input if an error occurs during processing.
        None: If the input is None or cannot be processed.
    """
    try:
        parts: Any = gnomad.split("&")
        if parts:
            return float(max(parts))
    except Exception as e:
        logging.debug(f"Error processing gnomad '{gnomad}': {e}")
        return gnomad


def pick_af_fields(var: dict) -> Dict[str, Any]:
    """
    Extracts allele frequency (AF) fields from a variant dictionary and returns
    them as a dictionary with specific keys.

    The function prioritizes allele frequency data in the following order:
    1. gnomAD_AF (if available, also includes gnomAD_MAX if present)
    2. gnomADg_AF (if gnomAD_AF is not available, also includes gnomAD_MAX if present)
    3. ExAC_MAF
    4. GMAF (Thousand Genomes)

    Args:
        var (dict): A dictionary representing a variant. It is expected to have
            the following structure:
            - "ALT": The alternate allele.
            - "INFO": A dictionary containing:
                - "CSQ": A list of dictionaries with allele frequency fields such as:
                    - "ExAC_MAF": ExAC allele frequency.
                    - "GMAF": Thousand Genomes allele frequency.
                    - "gnomAD_AF": gnomAD allele frequency.
                    - "gnomADg_AF": gnomAD genome allele frequency.
                    - "MAX_AF": Maximum allele frequency in gnomAD.

    Returns:
        Dict[str, Any]: A dictionary with the following keys:
            - "gnomad_frequency": The gnomAD allele frequency (or gnomAD genome frequency if gnomAD_AF is unavailable).
            - "gnomad_max": The maximum allele frequency in gnomAD (if available).
            - "exac_frequency": The ExAC allele frequency (if available).
            - "thousandG_frequency": The Thousand Genomes allele frequency (if available).
    """
    af: Dict[str, str] = {
        "gnomad_frequency": "",
        "gnomad_max": "",
        "exac_frequency": "",
        "thousandG_frequency": "",
    }
    allele: Any = var["ALT"]
    exac: float | Literal[0] = parse_allele_freq(
        var["INFO"]["CSQ"][0].get("ExAC_MAF"), allele
    )
    thousand_g: Any = parse_allele_freq(var["INFO"]["CSQ"][0].get("GMAF"), allele)
    gnomad: Any = var["INFO"]["CSQ"][0].get("gnomAD_AF", 0)
    gnomad_genome: Any = var["INFO"]["CSQ"][0].get("gnomADg_AF", 0)
    gnomad_max: Any = var["INFO"]["CSQ"][0].get("MAX_AF", 0)

    if gnomad:
        af["gnomad_frequency"] = max_gnomad(gnomad)
        if gnomad_max:
            af["gnomad_max"] = gnomad_max
    elif gnomad_genome:
        af["gnomad_frequency"] = gnomad_genome
        if gnomad_max:
            af["gnomad_max"] = gnomad_max
    if exac:
        af["exac_frequency"] = exac
    if thousand_g:
        af["thousandG_frequency"] = thousand_g
    return af


def read_mane(txt_gz: str) -> Dict[str, Dict[str, str]]:
    """
    Reads a gzipped tab-delimited file and creates a mapping of MANE (Matched Annotation
    from NCBI and EMBL-EBI) data to RefSeq and Ensembl identifiers.

    Args:
        txt_gz (str): Path to the gzipped tab-delimited file containing MANE data.

    Returns:
        Dict[str, Dict[str, str]]: A dictionary where each key is an Ensembl gene ID,
        and the value is another dictionary with:
            - "refseq": The RefSeq nucleotide ID (without version number).
            - "ensembl": The Ensembl nucleotide ID (without version number).

    The input file is expected to have the following columns:
        - "RefSeq_nuc": RefSeq nucleotide ID (e.g., "NM_001256789.1").
        - "Ensembl_nuc": Ensembl nucleotide ID (e.g., "ENST00000380152.8").
        - "Ensembl_Gene": Ensembl gene ID (e.g., "ENSG00000139618").
    """
    mane: Dict[str, Dict[str, str]] = {}
    with gzip.open(txt_gz, "rt") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for line in rdr:
            refseq: str | Any = line["RefSeq_nuc"].split(".")[0]
            ensembl: str | Any = line["Ensembl_nuc"].split(".")[0]
            gene: str | Any = line["Ensembl_Gene"].split(".")[0]
            mane[gene] = {"refseq": refseq, "ensembl": ensembl}
    return mane


def get_data_counts(
    preload: Dict[str, Any], sample_name: str, sid: str
) -> Dict[str, int | bool]:
    """
    Logs and returns the counts of various data types present in the preload dictionary.

    Args:
        preload (Dict[str, Any]): A dictionary containing preloaded data with keys
            such as 'snvs', 'cnvs', 'biomarkers', etc.
        sample_name (str): The name of the sample being processed.
        sid (str): The sample ID.

    Returns:
        Dict[str, int | bool]: A dictionary where each key corresponds to a data type
            in the preload dictionary, and the value is either the count of items
            (if the value is a list) or a boolean indicating presence (if the value
            is not a list).
    """
    data_counts: Dict[str, int | bool] = {
        key: (
            len(preload[key]) if isinstance(preload[key], list) else bool(preload[key])
        )
        for key in preload
    }
    logging.info(f"Data counts for sample {sample_name} ({sid}): {data_counts}")
    return data_counts


# --------------------------
# Parser Strategy
# --------------------------
class Parser(Protocol):
    """
    Parser interface for handling different data types.

    This interface defines a contract for parsers that process input arguments
    and return a dictionary containing preloaded data. Each parser implementation
    must provide its own logic for the `parse` method.
    """

    def parse(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """
        Processes the input arguments and returns a dictionary containing
        preloaded data. The returned dictionary may include any of the
        following keys:
            - snvs: Single nucleotide variants
            - cnvs: Copy number variants
            - biomarkers: Biomarker data
            - transloc: Translocations
            - lowcov: Low coverage regions
            - cov: Coverage data
            - fusions: Gene fusions
            - rna_expr: RNA expression data
            - rna_class: RNA classification data
            - qc: Quality control metrics
        """
        ...


@dataclass
class DnaParser:
    """
    DnaParser

    A class designed to parse and preprocess DNA data types, including VCF files, CNV JSON, Biomarkers JSON,
    DNA translocations, and coverage data. The parser processes and validates input data, extracts relevant
    information, and prepares it for downstream analysis.

    Attributes:
        canonical (Dict[str, str]): A dictionary mapping canonical gene/transcript information.

    Usage:
        This class is intended for use in pipelines that require preprocessing of DNA data for analysis.
        It ensures that the input data is validated, formatted, and enriched with relevant annotations
        before being passed to downstream tools or workflows.
    """

    canonical: Dict[str, str]

    def parse(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parses the input arguments and processes the DNA data files. Supports VCF, CNV, Biomarkers, Coverage,
        and translocations. Validates the presence of required files and ensures consistency in the input.

        Args:
            args (Dict[str, Any]): A dictionary containing input arguments and file paths for various data types.
                Supported keys include:
                    - vcf_files: Path to the VCF file containing SNVs.
                    - cnv: Path to the CNV JSON file.
                    - biomarkers: Path to the Biomarkers JSON file.
                    - transloc: Path to the DNA translocations VCF file.
                    - lowcov: Path to the low coverage BED file.
                    - cov: Path to the coverage JSON file.
                    - name: Sample name (optional).
                    - assay: Assay type (optional).
                    - update: Flag indicating whether to update existing data (optional).

        Returns:
            Dict[str, Any]: A dictionary containing preloaded data. The keys may include:
                - snvs: List of parsed SNVs (if VCF file is provided).
                - cnvs: List of parsed CNVs (if CNV JSON file is provided).
                - biomarkers: List of parsed biomarkers (if Biomarkers JSON file is provided).
                - transloc: List of parsed translocations (if translocations VCF file is provided
                - lowcov: List of parsed low coverage regions (if low coverage BED file is provided).
                - cov: Parsed coverage data (if coverage JSON file is provided).

        Raises:
            FileNotFoundError: If any of the required files (VCF, CNV JSON, Biomarkers JSON, translocations VCF, low coverage BED, coverage JSON) are missing.
            ValueError: If both lowcov and cov options are provided simultaneously.

        Notes:
            - The method uses helper functions to validate file existence and parse specific data types.
            - The method supports an "update" flag to skip parsing SNVs if the value is
                "no_update".
            - The method ensures that only one of lowcov or cov options is provided to avoid conflicts.
        """
        preload: Dict[str, Any] = {}
        vcf: Any | None = args.get("vcf_files")
        if not args.get("update") or not _is_no_update(vcf):
            _require_exists("VCF", vcf)
            preload["snvs"] = self._parse_snvs_only(vcf, args.get("assay"))

        if "cnv" in args:
            _require_exists("CNV JSON", args["cnv"])
            with open(args["cnv"], "r") as f:
                cnv_dict: Any = json.load(f)
            preload["cnvs"] = [cnv_dict[k] for k in cnv_dict]

        if "biomarkers" in args:
            _require_exists("Biomarkers JSON", args["biomarkers"])
            with open(args["biomarkers"], "r") as f:
                preload["biomarkers"] = json.load(f)

        if "transloc" in args:
            _require_exists("DNA translocations VCF", args["transloc"])
            preload["transloc"] = self._parse_transloc_only(args["transloc"])

        if "lowcov" in args and "cov" in args:
            raise ValueError("Both lowcov and cov present; choose one.")
        if "lowcov" in args:
            _require_exists("Low coverage BED", args["lowcov"])
            preload["lowcov"] = self._parse_lowcov_only(
                args["lowcov"], args.get("name")
            )
        if "cov" in args:
            _require_exists("Coverage JSON", args["cov"])
            with open(args["cov"], "r") as f:
                preload["cov"] = json.load(f)

        return preload

    def _parse_snvs_only(self, infile: str, assay: Optional[str]) -> list:
        """
        Parses single nucleotide variants (SNVs) from a VCF file. Filters variants based on
        specific criteria, processes variant annotations, and prepares the data for further analysis.
        Args:
            infile (str): Path to the VCF file containing SNVs.
            assay (Optional[str]): The assay type, which may influence filtering criteria.
        Returns:
            list: A list of dictionaries, each representing a parsed and processed variant.
        Notes:
            - The method uses the `pysam` library to read and iterate through the Variants
            - Variants are filtered based on the presence of specific annotations and criteria.
            - The method processes variant annotations, including consequence annotations,
            allele frequencies, and other relevant information.
            - The method ensures that the variant data is formatted and enriched with
            relevant annotations before being returned.
        """
        filtered_data: list = []
        vcf_object = VariantFile(infile)
        for var in vcf_object.fetch():
            var_dict: dict[str, Any] = cmdvcf.parse_variant(var, vcf_object.header)
            var_csq = var_dict["INFO"]["CSQ"]
            # removes the variant all together if the csq only has experimental transcripts but retains if the any of those experimental transcripts are in the list of genes
            if var_csq:
                all_features: list = [c.get("Feature") for c in var_csq]
                all_X_genes: list = [
                    c.get("SYMBOL")
                    for c in var_csq
                    if c.get("Feature", "").startswith("X")
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
            # find floats in CSQ, make sure they are saved into coyote as such
            var_dict = emulate_perl(var_dict)
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
                slim_csq, self.canonical
            )  # How slow is this?

            var_dict["INFO"]["CSQ"] = selected_transcript_removal(
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
            var_dict["INFO"]["variant_callers"] = var_dict["INFO"][
                "variant_callers"
            ].split("|")
            var_dict["FILTER"] = var_dict["FILTER"].split(";")

            # Keep parity with the legacy Perl importer:
            # skip variants that failed NVAF/LONGDEL or failed any PON filter.
            filters = set(var_dict["FILTER"])
            if "FAIL_NVAF" in filters or "FAIL_LONGDEL" in filters:
                continue
            if any(f.startswith("FAIL_PON_") for f in filters):
                continue
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
                    exit(
                        "not a valid VCF, should be aggregated by AF(VAF), VD AD and GT"
                    )
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
        return filtered_data

    def _parse_lowcov_only(self, lowcov_bed: str, sample_name: Optional[str]) -> list:
        """
        Parses low coverage regions from a BED file. Extracts and formats information about genomic regions
        with low coverage, including average coverage and amplicon details.
        Args:
            lowcov_bed (str): Path to the BED file containing low coverage regions.
            sample_name (Optional[str]): The name of the sample associated with the low coverage data.
        Returns:
            list: A list of dictionaries, each representing a low coverage region with relevant details.
        Notes:
            - The method reads the BED file and extracts information such as chromosome, start and end positions
            average coverage, and amplicon details.
            - The sample name is added to each low coverage region entry for identification.
            - The method ensures that the start and end positions are stored as integers and the average coverage
            is stored as a float.
        """
        lowcov_data: list = []
        with open(lowcov_bed, "r") as f:
            lowcov_dict = csv.DictReader(
                f,
                delimiter="\t",
                fieldnames=["chr", "start", "end", "avg_cov", "amplicon"],
            )
            for row in lowcov_dict:
                row["sample"] = sample_name
                row["start"] = int(row["start"])
                row["end"] = int(row["end"])
                row["avg_cov"] = float(row["avg_cov"])
                lowcov_data.append(row)
        return lowcov_data

    def _parse_transloc_only(self, infile: str) -> list:
        """
        Parses DNA translocations from a VCF (Variant Call Format) file, filtering and processing
        translocation annotations such as gene fusions and bidirectional fusions. Additionally,
        integrates MANE (Matched Annotation from NCBI and EMBL-EBI) annotations where applicable.

        Args:
            infile (str): Path to the input VCF file containing variant data.

        Returns:
            list: A list of dictionaries representing filtered and processed variant data. Each
            dictionary contains variant information, including annotations and MANE data if available.

        Functionality:
            - Reads MANE annotations from a configuration file.
            - Filters out variants with "<" in their ALT field (e.g., deletions and duplications).
            - Processes variant annotations to identify gene fusions and bidirectional fusions.
            - Matches MANE annotations for both genes in a pair and adds them to the variant data.
            - Removes dot notation from SNPeff annotations for database compatibility.
            - Constructs a new annotation structure for each variant and appends it to the output
            if it meets the filtering criteria.
        """
        mane: Dict[str, Dict[str, str]] = read_mane(config.mane)
        filtered_data: list = []
        vcf_object = VariantFile(infile)
        fusion_buffer: Dict[str, Any] = {}

        for var in vcf_object.fetch():
            var_dict: dict[str, Any] = cmdvcf.parse_variant(var, vcf_object.header)

            keep_variant = 0
            mane_select = {}
            all_new_ann = []
            add_mane = 0
            svtype = var.info.get("SVTYPE", "None")
            
            if svtype not in {"BND", "DUP", "DEL"}:
                continue
            
            if var.filter.keys() != ['PASS']:
                continue
            
            for ann in var_dict["INFO"]["ANN"]:
                if any(anno in ["gene_fusion", "bidirectional_gene_fusion", "frameshift_variant","feature_fusion"] for anno in ann.get("Annotation", [])):

                    if ann.get("Annotation") == ['feature_fusion']:
                        if ann.get("Feature_Type") == "CUSTOM&sorted":
                            feature_genes, feature_id = ann["Feature_ID"].split("_", 1)
                            #print (feature_genes, feature_id)                
                            var_id = var.id
                            mate_id = var.info.get("MATEID")
                            hgvs = ann.get("HGVS.c")

                            if isinstance(mate_id, (tuple, list)):
                                mate_id = mate_id[0] 

                            fusion_buffer[var_id] = {
                                "gene": feature_genes,
                                "id": feature_id,
                                "ann": ann,
                                "hgvs": hgvs,
                                "mate": mate_id
                                }

                            #print (mate_id)

                            if mate_id in fusion_buffer:
                                #print(fusion_buffer)
                                f1 = fusion_buffer[mate_id]
                                f2 = {
                                    "gene": feature_genes,
                                    "id": feature_id,
                                    "ann": ann,
                                    "hgvs": hgvs,
                                }

                                halves = sorted([f1, f2], key=lambda x: x["gene"].upper())

                                g1, g2 = halves[0]["gene"], halves[1]["gene"]
                                id1, id2 = halves[0]["id"], halves[1]["id"]
                                hgvs1,hgvs2 = halves[0]["hgvs"], halves[1]["hgvs"]

                        
                                combined_ann = {
                                    "Allele": var.alleles[0],
                                    "Annotation": ["feature_fusion"],
                                    "Annotation_Impact": "LOW",
                                    "Gene_Name": f"{g1}&{g2}",
                                    "Gene_ID": f"{id1}&{id2}",
                                    "Feature_Type": "CUSTOM&sorted",
                                    "Feature_ID": f"{g1}_{id1}&{g2}_{id2}",
                                    "Transcript_BioType": "",
                                    "Rank": "",
                                    "HGVS.c": f"{hgvs1};{hgvs2}",   # optional symbolic fusion notation
                                    "HGVS.p": "",
                                    "cDNA.pos / cDNA.length": "",
                                    "CDS.pos / CDS.length": "",
                                    "AA.pos / AA.length": "",
                                    "Distance": "",
                                    "ERRORS / WARNINGS / INFO": "",
                                }

                                #print(combined_ann)
                                ann = combined_ann
                                genes = ann["Gene_ID"].split("&")
                                #print(ann.get("Annotation"), ann.get("Feature_Type"), ann.get("Gene_Name"), genes)

                                del fusion_buffer[mate_id]
                                del fusion_buffer[var_id]
                            else:
                                continue
                        else:
                            continue

                    else:
                        genes = ann["Gene_ID"].split("&")
                        #print(ann.get("Annotation"), ann.get("Feature_Type"), ann.get("Gene_Name"), genes)

                    keep_variant = 1
                    n_mane = 0

                    if len(genes) > 2:
                        continue

                    for gene in genes:
                        enst = mane.get(gene, {}).get("ensembl", "NO_MANE_TRANSCRIPT")
                        if enst in ann.get("HGVS.p", ""):
                            n_mane += 1

                    new_ann = {key.replace(".", ""): ann[key] for key in ann}
                    all_new_ann.append(new_ann)
                
                    if n_mane > 0 and n_mane == len(genes) and isinstance(mane_select, dict) and not mane_select:
                        mane_select = new_ann
                        add_mane = 1
                    
            #print (mane_select)
            del var_dict["INFO"]["ANN"]
            var_dict["INFO"]["ANN"] = all_new_ann
            
            if add_mane:
                var_dict["INFO"]["MANE_ANN"] = mane_select
            if keep_variant:
                filtered_data.append(var_dict)
                
        return filtered_data


@dataclass
class RnaParser:
    """
    A class designed to parse and preprocess RNA data types, including gene fusions,
    RNA expression data, RNA classification data, and quality control metrics. The parser
    processes and validates input data, extracts relevant information, and prepares it for
    downstream analysis.

    Usage:
        This class is intended for use in pipelines that require preprocessing of RNA data
        for analysis. It ensures that the input data is validated, formatted, and enriched
        with relevant annotations before being passed to downstream tools or workflows.
    """

    def parse(self, args: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parses the input arguments and processes the RNA data files. Supports gene fusions,
        RNA expression data, RNA classification data, and quality control metrics. Validates
        the presence of required files and ensures consistency in the input.

        Args:
            args (Dict[str, Any]): A dictionary containing input arguments and file paths
            for various data types. Supported keys include:
                - fusion_files: Path to the gene fusions JSON file.
                - expression_path: Path to the RNA expression JSON file.
                - classification_path: Path to the RNA classification JSON file.
                - qc: Path to the quality control JSON file.
                - name: Sample name (optional).
                - update: Flag indicating whether to update existing data (optional).

        Returns:
            Dict[str, Any]: A dictionary containing preloaded data. The keys may include:
                - fusions: List of parsed gene fusions (if fusion JSON file is provided).
                - rna_expr: RNA expression data (if RNA expression JSON file is provided).
                - rna_class: RNA classification data (if RNA classification JSON file is provided).
                - qc: Quality control metrics (if QC JSON file is provided).

        Raises:
            FileNotFoundError: If any of the required files (gene fusions JSON, RNA
            expression JSON, RNA classification JSON, QC JSON) are missing.

        Notes:
            - The method uses helper functions to validate file existence and parse
            specific data types.
            - The method supports an "update" flag to skip parsing fusions if the value
                is "no_update".
        """
        preload: Dict[str, Any] = {}
        fusions: Any | None = args.get("fusion_files")
        if not args.get("update") or not _is_no_update(fusions):
            _require_exists("Fusions JSON", fusions)
            with open(fusions, "r") as f:
                preload["fusions"] = json.load(f)

        if "expression_path" in args:
            _require_exists("Expression JSON", args["expression_path"])
            with open(args["expression_path"], "r") as f:
                preload["rna_expr"] = json.load(f)
        if "classification_path" in args:
            _require_exists("Classification JSON", args["classification_path"])
            with open(args["classification_path"], "r") as f:
                preload["rna_class"] = json.load(f)
        if "qc" in args:
            _require_exists("QC JSON", args["qc"])
            with open(args["qc"], "r") as f:
                preload["rna_qc"] = json.load(f)

        return preload


# --------------------------
# Writer (single write path)
# --------------------------
class DependentWriter:
    """
    A class responsible for writing preloaded data to a database, handling various
    collections based on the type of data. The writer ensures that existing entries
    for a sample can be deleted if specified, and it manages the insertion of new
    data into the appropriate collections.
    """

    def __init__(self, repos: Repos) -> None:
        self.repos: Repos = repos
        self.collection_binds: Dict[str, str] = {
            # preload key : collection name
            "snvs": "variants",
            "cnvs": "cnvs",
            "biomarkers": "biomarkers",
            "transloc": "transloc",
            "lowcov": "coverage",
            "cov": "panel_cov",
            "fusions": "fusions",
            "rna_expr": "rna_expression",
            "rna_class": "rna_classification",
            "rna_qc": "rna_qc",
        }

    def write(
        self,
        preload: Dict[str, Any],
        sample_id: ObjectId,
        sample_name: str,
        delete_existing: bool = False,
    ) -> None:
        """
        Writes preloaded data to the database, handling various collections based on
        the type of data. Existing entries for a sample can be deleted if specified.

        Args:
            preload (Dict[str, Any]): A dictionary containing preloaded data with keys
                such as 'snvs', 'cnvs', 'biomarkers', etc.
            sample_id (ObjectId): The unique identifier for the sample.
            sample_name (str): The name of the sample being processed.
            delete_existing (bool): Flag indicating whether to delete existing entries
                for the sample before inserting new data. Defaults to False.

        Raises:
            TypeError: If the data types of the payloads do not match the expected types
                (e.g., expecting a list of dictionaries but receiving a different type).
        """
        sid = str(sample_id)

        def wipe(col: str) -> None:
            """
            Deletes existing entries for the sample in the specified collection.
            Args:
                col (str): The name of the collection from which to delete entries.
            """
            if delete_existing:
                logging.info(
                    f"Wiping existing entries in {col} for sample {sample_name}({sid})"
                )
                self.repos.col(col).delete_many({"SAMPLE_ID": sid})

        # Helpful debug to catch type mismatches early
        logging.debug(
            "preload types: %s", {k: type(v).__name__ for k, v in preload.items()}
        )

        for key, col in self.collection_binds.items():
            if key not in preload:
                continue

            payload: Any = preload[key]
            wipe(col)

            # dict-shaped payloads
            if key in ("biomarkers", "cov", "rna_expr", "rna_class", "rna_qc"):
                if not isinstance(payload, dict):
                    raise TypeError(
                        f"{key} expected dict, got {type(payload).__name__}"
                    )
                payload["SAMPLE_ID"] = sid
                if key == "cov":
                    payload["sample"] = sample_name
                self.repos.col(col).insert_one(payload)
                logging.info(f"Inserted {key} data for sample {sample_name} ({sid})")
                continue

            # list/tuple-shaped payloads
            if not isinstance(payload, (list, tuple)):
                raise TypeError(
                    f"{key} expected list of dicts, got {type(payload).__name__}"
                )
            for rec in payload:
                if not isinstance(rec, dict):
                    raise TypeError(
                        f"{key} contained {type(rec).__name__}, expected dict"
                    )
                rec["SAMPLE_ID"] = sid
            if payload:
                self.repos.col(col).insert_many(payload)
                logging.info(
                    f"Inserted {len(payload)} {key} records for sample {sample_name} ({sid})"
                )

        # Embedded blobs on the sample doc

        # if "rna_expr" in preload:
        #     self.repos.samples.update_one(
        #         {"_id": sid}, {"$set": {"expr": preload["rna_expr"]}}
        #     )
        #     logging.info(
        #         f"Inserted RNA expression data for sample {sample_name} ({sid})"
        #     )
        # if "rna_class" in preload:
        #     self.repos.samples.update_one(
        #         {"_id": sid}, {"$set": {"classification": preload["rna_class"]}}
        #     )
        #     logging.info(
        #         f"Inserted RNA classification data for sample {sample_name} ({sid})"
        #     )
        # if "rna_qc" in preload:
        #     self.repos.samples.update_one(
        #         {"_id": sid}, {"$set": {"QC": [preload["rna_qc"]]}}
        #     )
        #     logging.info(f"Inserted QC data for sample {sample_name} ({sid})")


# --------------------------
# Ingestion Orchestrator
# --------------------------
@contextmanager
def cleanup_on_error(repos: Repos, sample_id: ObjectId) -> Generator[None, Any, None]:
    """
    This context manager is designed to handle errors that may occur during the ingestion
    process. If an exception is raised within the context, it ensures that all related
    database entries associated with the given `sample_id` are cleaned up to maintain
    data consistency.

    Args:
        repos (Repos): An instance of the `Repos` class, which provides access to the
            database collections.
        sample_id (ObjectId): The unique identifier of the sample being ingested.

    Yields:
        None: The context manager does not yield any specific value.

    Raises:
        Exception: Re-raises the original exception after performing cleanup operations.

    Behavior:
        - Logs an error message indicating the failure of ingestion and the initiation
        of cleanup.
        - Iterates through a predefined list of database collections (`variants`, `cnvs`,
        `biomarkers`, `transloc`, `coverage`, `panel_cov`, `fusions`) and deletes all
        entries associated with the `sample_id`.
        - Deletes the sample entry from the `samples` collection.
        - Suppresses any exceptions that occur during the cleanup process to ensure
        that the original exception is re-raised.
    """
    try:
        yield
    except Exception as e:
        sid = str(sample_id)
        logging.error(f"Ingest failed, cleaning up: {e}")
        for col in (
            "variants",
            "cnvs",
            "biomarkers",
            "transloc",
            "coverage",
            "panel_cov",
            "fusions",
            "rna_expression",
            "rna_classification",
            "rna_qc",
        ):
            try:
                repos.col(col).delete_many({"SAMPLE_ID": sid})
                logging.info(f"Successfully cleaned up {col} for sample {sid}.")
            except Exception as e:
                logging.error(f"Failed to clean up {col} for {sid}: {e}")
        try:
            repos.samples.delete_one({"_id": sample_id})
            logging.info(
                f"Successfully cleaned up sample {sid} from samples collection."
            )
        except Exception as es:
            logging.error(f"Failed to clean up samples for {sid}: {es}")
        finally:
            logging.error(f"Raising exception after cleanup for sample {sid}.")
        raise


class IngestionService:
    """
    A service class responsible for managing the ingestion and updating of sample data
    into a repository. This class handles the creation of new samples, updating existing
    samples, and ensuring data consistency between DNA and RNA samples.

    Attributes:
        repos (Repos): The repository interface for interacting with the database.
        parser (Parser): A parser instance for processing input arguments.
        writer (DependentWriter): A writer instance for managing dependent data.
        canonical (Dict[str, str]): An optional dictionary for canonical mappings.
    """

    def __init__(
        self, repos: Repos, parser: Parser, canonical: Optional[Dict[str, str]] = None
    ) -> None:
        self.repos: Repos = repos
        self.parser: Parser = parser
        self.writer = DependentWriter(repos)
        self.canonical: Dict[str, str] = canonical or {}

    def ingest_new(self, args: Dict[str, Any]) -> ObjectId:
        """
        This method performs the following steps:
        1. Parses the input arguments to prepare the data for ingestion.
        2. Generates a unique name for the sample, optionally incrementing if required.
        3. Writes dependent data to the database using a pre-generated sample ID.
        4. Constructs metadata for the sample, including data counts and timestamps.
        5. Inserts the sample metadata into the repository.

        Args:
            args (Dict[str, Any]): A dictionary of arguments containing the data to be ingested.
                Expected keys include:
                - "name" (str): The base name for the sample.
                - "increment" (Optional[bool]): Whether to increment the name if it already exists.

        Returns:
            ObjectId: The unique identifier of the newly ingested sample.

        Raises:
            Exception: If an error occurs during the ingestion process, the operation is rolled back.

        Notes:
            - The `cleanup_on_error` context manager ensures that any partial changes are reverted
            in case of an error during the ingestion process.
            - The `data_counts` dictionary provides a summary of the loaded data for the sample.
        """
        preload: Dict[str, Any] = self.parser.parse(args)
        name: str = what_id(args["name"], args.get("increment"), self.repos.samples)
        sample_id = ObjectId()
        with cleanup_on_error(self.repos, sample_id):
            # dependents first (use pre-generated _id)
            self.writer.write(preload, sample_id, name, delete_existing=False)

            # sample last
            meta: Dict[str, Any] = build_sample_meta_dict(args)

            # Add loaded data counts to sample
            data_counts: Dict[str, int | bool] = get_data_counts(
                preload, name, str(sample_id)
            )

            meta.update(
                {
                    "_id": sample_id,
                    "name": name,
                    "data_counts": data_counts,
                    "time_added": datetime.utcnow(),
                    "ingest_status": "ready",
                }
            )
            self.repos.samples.insert_one(meta)
            logging.info(f"Ingest complete for {name} ({sample_id}).")
        return sample_id

    def update_existing(self, args: Dict[str, Any]) -> ObjectId:
        """
        Updates an existing sample in the database with new data provided in the `args` dictionary.

        Args:
            args (Dict[str, Any]): A dictionary containing the data to update the sample with.
                Must include the "name" key to identify the sample. Additional keys may include
                data files and metadata.

        Returns:
            ObjectId: The unique identifier of the updated sample in the database.

        Raises:
            SystemExit: If the sample cannot be found in the database, if there is a mismatch
                between the existing sample type (DNA/RNA) and the provided data, or if the
                sample type cannot be determined.

        Workflow:
            1. Checks if the sample exists in the database using the "name" key in `args`.
            2. Ensures consistency between the existing sample type (DNA/RNA) and the provided data.
            3. Parses the input data for validation and preloading.
            4. Updates the sample's metadata, blocking changes to specific fields like "assay".
            5. Calculates updated data counts for the sample.
            6. Updates the sample's ingest status and data counts in the database.
            7. Deletes and rewrites dependent data associated with the sample.
            8. Logs the completion of the update process.

        Notes:
            - The function enforces strict DNA/RNA consistency and allows bypassing updates
                for mandatory files using the "no_update" keyword.
            - The `data_typer` function is used to determine the type of data (DNA/RNA) in `args`.
            - The `parser`, `_update_meta`, `get_data_counts`, and `writer` are assumed to be
                components of the class or external utilities used for parsing, updating, and writing data.
        """
        doc = self.repos.samples.find_one({"name": args["name"]})
        if not doc:
            exit("Cannot find case in database, will not update anything. Bye!")

        sample_id: ObjectId = doc["_id"]

        # Enforce DNA/RNA consistency and allow "no_update" bypass for mandatory files
        if "fusion_files" in doc:  # existing RNA
            if data_typer(args) == "DNA":
                exit("you are trying to add DNA data to a RNA sample. BAD PERSON!")
            if "fusion_files" not in args:
                args["fusion_files"] = "no_update"
        elif "vcf_files" in doc:  # existing DNA
            if data_typer(args) == "RNA":
                exit("you are trying to add RNA data to a DNA sample. BAD PERSON!")
            if "vcf_files" not in args:
                args["vcf_files"] = "no_update"
        else:
            exit("update function could not determine if the case is DNA or RNA")

        # Parse first (preflight)
        preload: Dict[str, Any] = self.parser.parse(args)

        # Update metadata (block assay change)
        self._update_meta(
            build_sample_meta_dict(args), sample_id, block_fields={"assay"}
        )

        # Get updated data counts
        data_counts: Dict[str, int | bool] = get_data_counts(
            preload, args["name"], str(sample_id)
        )
        # Ensure ready
        self.repos.samples.update_one(
            {"_id": sample_id},
            {"$set": {"ingest_status": "ready", "data_counts": data_counts}},
        )

        # Rewrite dependents (delete + insert)
        self.writer.write(preload, sample_id, args["name"], delete_existing=True)

        logging.info(f"Update complete for {args['name']} ({sample_id}).")
        return sample_id

    def _update_meta(
        self, meta: Dict[str, Any], sample_id: ObjectId, block_fields: set[str]
    ) -> None:
        """
        Updates the metadata of a sample in the database.

        This method compares the provided metadata (`meta`) with the current metadata
        of the sample identified by `sample_id`. If there are differences, it updates
        the database accordingly. Certain fields specified in `block_fields` cannot
        be updated, and attempting to do so will result in an exit with an error message.

        Args:
            meta (Dict[str, Any]): A dictionary containing the metadata to update.
            sample_id (ObjectId): The unique identifier of the sample to update.
            block_fields (set[str]): A set of field names that are restricted from being updated.

        Returns:
            None

        Raises:
            SystemExit: If an attempt is made to update a field listed in `block_fields`.
        """
        current = self.repos.samples.find_one({"_id": sample_id}) or {}
        for k, v in meta.items():
            if k in current and current[k] != v:
                if k in block_fields:
                    exit(f"No support to update {k} as of yet")
                self.repos.samples.update_one({"_id": sample_id}, {"$set": {k: v}})
                logging.debug(f"changing {k}: from {current[k]} to {v}")
            elif k not in current:
                self.repos.samples.update_one({"_id": sample_id}, {"$set": {k: v}})
                logging.debug(f"adding {k}: {v} to sample")


# --------------------------
# Main entry
# --------------------------
def main(args) -> None:
    command: Any = args.command_selection
    if command == "load":
        args_dict: dict = {k: v for k, v in vars(args).items() if v is not None}
    elif command == "yaml":
        args_dict: Dict[str, Any] = validate_yaml(args.yaml_file)
        args_dict["update"] = args.update
        args_dict["increment"] = args.increment
        args_dict["dev"] = getattr(args, "dev", False)
        args_dict["debug_logger"] = getattr(args, "debug_logger", False)
    else:
        exit("Unknown command")

    # Decide DNA vs RNA
    dtype: str | None = data_typer(args_dict)
    if dtype not in ("DNA", "RNA"):
        exit("Could not determine data type (DNA/RNA) from inputs.")

    # DB + canonical
    repos: Repos = Repos.from_args(args_dict)
    canonical: dict = {}
    if dtype == "DNA":
        # load canonical map once for this DB
        canonical: dict = {}
        for c in repos.col("refseq_canonical").find({}):
            canonical[c["gene"]] = c["canonical"]

    # Parser strategy
    parser: Parser = DnaParser(canonical) if dtype == "DNA" else RnaParser()
    service = IngestionService(repos=repos, parser=parser, canonical=canonical)

    if args_dict.get("update"):
        service.update_existing(args_dict)
    else:
        service.ingest_new(args_dict)


if __name__ == "__main__":
    parser: ArgumentParser = cli_parser()
    args: Namespace = parser.parse_args()
    if not args.command_selection:
        parser.print_help()
        sys.exit(2)
    setup_logging(debug=args.debug_logger)
    main(args)
