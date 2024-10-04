assay = {
    'solid_GMSv3': {
        'filters': {
            'gnomAD_AF': 0.05,
            'FAIL_PON': 1,
            'FAIL_NVAF': 1,
            'FAIL_LONGDEL': 1
        }
    },
    'myeloid_GMSv1': {
        'filters': {
            'gnomAD_AF': 0.05,
            'FAIL_PON': 1,
            'FAIL_NVAF': 1,
            'FAIL_LONGDEL': 1
        }
    },
    'tumwgs': {
        'filters': {
            'gnomAD_AF': 0.05,
            'FAIL_PON': 1,
            'FAIL_NVAF': 1,
            'FAIL_LONGDEL': 1
        }
    },
    'PARP_inhib': {
        'filters': {
            'gnomAD_AF': 0.05,
            'FAIL_PON': 1,
            'FAIL_NVAF': 1,
            'FAIL_LONGDEL': 1
        }
    },
    'default': {
        'filters': {
            'gnomAD_AF': 0.05,
            'FAIL_PON': 1,
            'FAIL_NVAF': 1,
            'FAIL_LONGDEL': 1
        }
    },
}
coyote_csq = {
    "AA_AF": 1,
    "AFR_AF": 1,
    "Amino_acids": 1,
    "AMR_AF": 1,
    "APPRIS": 1,
    "BAM_EDIT": 1,
    "CADD_RAW": 1,
    "CCDS": 1,
    "cDNA_position": 1,
    "CDS_position": 1,
    "Codons": 1,
    "DISTANCE": 1,
    "DOMAINS": 1,
    "EA_AF": 1,
    "EAS_AF": 1,
    "ENSP": 1,
    "EUR_AF": 1,
    "Feature_type": 1,
    "GENE_PHENO": 1,
    "GIVEN_REF": 1,
    "gnomADg": 1,
    "gnomADg_AF_popmax": 1,
    "gnomADg_popmax": 1,
    "gnomAD_OTH_AF": 1,
    "HGVS_OFFSET": 1,
    "HIGH_INF_POS": 1,
    "LoFtool": 1,
    "MAX_AF": 1,
    "MAX_AF_POPS": 1,
    "miRNA": 1,
    "MOTIF_NAME": 1,
    "MOTIF_POS": 1,
    "MOTIF_SCORE_CHANGE": 1,
    "PHENO": 1,
    "Protein_position": 1,
    "PUBMED": 1,
    "REFSEQ_MATCH": 1,
    "SAS_AF": 1,
    "SOMATIC": 1,
    "STRAND": 1,
    "SWISSPROT": 1,
    "SYMBOL_SOURCE": 1,
    "TREMBL": 1,
    "UNIPARC": 1,
    "USED_REF": 1
}
mane = "/data/bnf/dev/viktor/cmdvcf/resources/MANE.GRCh38.v0.9.summary.txt.gz"
mongo = {
    "uri": "mtlucmds1.lund.skane.se",
    "dbname": "coyote_dev",
}

data_types = {
    'expression_path' : "RNA",
    'classification_path' : "RNA",
    'fusion_files' : "RNA",
    'qc' : "RNA", 
    'biomarkers' : "DNA",
    'vcf_files' : "DNA",
    'lowcov' : "DNA",
    'transloc' : "DNA",
    'cnvprofile' : "DNA",
    'cnv' : "DNA",
}