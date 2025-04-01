assay = {
    "solid_GMSv3": {
        "filters": {"gnomAD_AF": 0.05, "FAIL_PON": 1, "FAIL_NVAF": 1, "FAIL_LONGDEL": 1}
    },
    "myeloid_GMSv1": {
        "filters": {"gnomAD_AF": 0.05, "FAIL_PON": 1, "FAIL_NVAF": 1, "FAIL_LONGDEL": 1}
    },
    "tumwgs": {
        "filters": {"gnomAD_AF": 0.05, "FAIL_PON": 1, "FAIL_NVAF": 1, "FAIL_LONGDEL": 1}
    },
    "PARP_inhib": {
        "filters": {"gnomAD_AF": 0.05, "FAIL_PON": 1, "FAIL_NVAF": 1, "FAIL_LONGDEL": 1}
    },
    "default": {
        "filters": {"gnomAD_AF": 0.05, "FAIL_PON": 1, "FAIL_NVAF": 1, "FAIL_LONGDEL": 1}
    },
    'GMSHem': {
        'filters': {'gnomAD_AF': 0.05,'FAIL_PON': 1,'FAIL_NVAF': 1,'FAIL_LONGDEL': 1}
    },
    'default': {
        'filters': {'gnomAD_AF': 0.05,'FAIL_PON': 1,'FAIL_NVAF': 1,'FAIL_LONGDEL': 1}
    },
}

mane = "/data/bnf/dev/viktor/cmdvcf/resources/MANE.GRCh38.v0.9.summary.txt.gz"
mongo = {
    "uri": "mtlucmds1.lund.skane.se",
    "dbname": "coyote_dev_3",
}

data_types = {
    "expression_path": "RNA",
    "classification_path": "RNA",
    "fusion_files": "RNA",
    "qc": "RNA",
    "biomarkers": "DNA",
    "vcf_files": "DNA",
    "lowcov": "DNA",
    "transloc": "DNA",
    "cnvprofile": "DNA",
    "cnv": "DNA",
}
