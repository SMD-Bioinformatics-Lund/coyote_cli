# Coyote input data specification

This document describes the input file shapes currently expected by `import_coyote_sample.py` when loading SNVs/indels, CNVs, biomarkers, and RNA fusions into Coyote.

The loader does not perform schema validation for most JSON inputs. Fields described as required are required by the code path itself or by the downstream Coyote data model inferred from the example files.

## Overview

| Data type | CLI/YAML key | File format | Loader target collection |
| --- | --- | --- | --- |
| SNVs/indels | `vcf_files` / `--vcf` | VCF, VEP annotated | `variants` |
| CNVs | `cnv` / `--cnv` | JSON object | `cnvs` |
| Biomarkers | `biomarkers` / `--biomarkers` | JSON object | `biomarkers` |
| RNA fusions | `fusion_files` / `--fusions` | JSON array | `fusions` |

For YAML loading, use the same keys shown above. See `resources/example.yaml` for a complete sample metadata example.

## SNVs and indels VCF

SNV/indel data must be supplied as a valid VCF readable by `pysam.VariantFile`. Coyote expects a VEP-annotated VCF with an INFO `CSQ` field and one or more sample genotype columns.

### Required record organization

Each variant record must contain:

| VCF field | Requirement |
| --- | --- |
| `CHROM` | Chromosome name. Stored as `CHROM`. |
| `POS` | 1-based position. Stored as `POS`. |
| `ID` | May be `.`. Stored as `ID`. |
| `REF` | Reference allele. Stored as `REF`. |
| `ALT` | Alternate allele. Stored as `ALT`. Prefer one ALT allele per record. |
| `QUAL` | Stored as `QUAL`; may be missing if valid VCF. |
| `FILTER` | Stored and split on `;` into a list. |
| `INFO/CSQ` | Required VEP transcript annotations. At least one transcript annotation must exist. |
| `INFO/variant_callers` | Required. Must be a string with callers separated by `|`, for example `vardict|freebayes`. |
| `FORMAT` | Must include `GT`, `VAF`, `VD`, and `DP`. |
| Sample columns | First sample column is treated as tumor/case. Any later sample columns are treated as control/normal. |

The sample column order is critical. The importer preserves VCF header sample order and assigns:

| Sample index | Coyote type |
| --- | --- |
| First sample column | `case` |
| Second and later sample columns | `control` |

For paired tumor/normal VCFs, the header must therefore end with tumor first and normal second:

```text
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT TumorSample NormalSample
```

For tumor-only VCFs, include one sample column. That single sample is loaded as the `case`.

### Required genotype FORMAT fields

The loader expects the following fields in every sample genotype:

| FORMAT field | Meaning | Example |
| --- | --- | --- |
| `GT` | Genotype. Pysam tuple values are converted to strings like `0/1`. | `0/1` |
| `VAF` | Alternate allele fraction. This is renamed to `AF` during import. | `0.5727` |
| `VD` | Alternate allele read count. | `1529` |
| `DP` | Total read depth. | `2670` |

Example FORMAT and sample values:

```text
GT:VAF:VD:DP    0/1:0.5727:1529:2670    0/1:0.01:2:2670
```

Although an older error message mentions `AF(VAF)` and `AD`, the current code uses `VAF`, `VD`, `DP`, and `GT`. Supplying only `AF` without `VAF` will fail later because the importer renames `VAF` to `AF`.

### Required INFO fields

The importer directly reads these INFO fields:

| INFO field | Requirement |
| --- | --- |
| `CSQ` | Required. Parsed according to the VEP `##INFO=<ID=CSQ...Format: ...>` header. |
| `variant_callers` | Required. Split on `|` and stored as a list. |
| `SVTYPE` | Optional. If present, copied to `INFO/TYPE`. |

### Required VEP CSQ header

The VCF header must define `CSQ` with a VEP-style `Format:` string. Coyote uses that header to map each pipe-delimited CSQ value into named transcript fields.

Minimum fields required for the current import logic:

| CSQ field | Why Coyote needs it |
| --- | --- |
| `Consequence` | Stored as a list split on `&`. |
| `IMPACT` | Used to select the displayed transcript by severity order: `HIGH`, `MODERATE`, `LOW`, `MODIFIER`. |
| `SYMBOL` | Used for gene summaries and canonical transcript lookup. |
| `Feature` | Transcript identifier. Used for transcript summaries and selected transcript tracking. |
| `BIOTYPE` | Used as fallback transcript selection, preferring `protein_coding`. |
| `CANONICAL` | Used as fallback transcript selection when `CANONICAL=YES`. |

Recommended fields for full Coyote display and searching:

| CSQ field | Use |
| --- | --- |
| `HGNC_ID` | Stored in transcript summaries. |
| `PolyPhen` | Stored in transcript summaries. |
| `SIFT` | Stored in transcript summaries. |
| `ENSP` | Stored in transcript summaries. |
| `INTRON` | Stored in transcript summaries. |
| `EXON` | Stored in transcript summaries. |
| `MANE_SELECT` | Stored as `MANE`. |
| `STRAND` | Stored in transcript summaries. |
| `CADD_PHRED` | Stored in transcript summaries. |
| `CLIN_SIG` | Stored in transcript summaries. |
| `VARIANT_CLASS` | Stored per transcript and as top-level `variant_class` from the first CSQ annotation. |
| `HGVSc` | Summarized into top-level `HGVSc`. The transcript prefix before `:` is removed. |
| `HGVSp` | Summarized into top-level `HGVSp`. The transcript/protein prefix before `:` is removed. |
| `Existing_variation` | Used to extract the first `rs...` dbSNP identifier. |
| `COSMIC` | Split on `&` into top-level `cosmic_ids`. |
| `PUBMED` | Split on `&` into top-level `pubmed_ids`. |
| `gnomAD_AF` | Used for top-level `gnomad_frequency`. Multi-values separated by `&` are reduced to the max. |
| `gnomADg_AF` | Fallback for top-level `gnomad_frequency`. |
| `MAX_AF` | Stored as top-level `gnomad_max` when gnomAD frequency is present. |
| `ExAC_MAF` | Used for top-level `exac_frequency` if formatted as `allele:frequency`. |
| `GMAF` | Used for top-level `thousandG_frequency` if formatted as `allele:frequency`. |
| `{d,gi,lu,cns,mm,co}hotspot_OID` | Non-empty values are collected into top-level `hotspots`. |

The example VCF in `resources/test.vcf` contains the expected broad CSQ shape.

### Transcript selection

Coyote keeps all slimmed transcript annotations in `INFO/CSQ` and also stores one selected transcript in `INFO/selected_CSQ`.

Selection order:

1. Iterate by `IMPACT` severity: `HIGH`, then `MODERATE`, then `LOW`, then `MODIFIER`.
2. Prefer a transcript matching the `refseq_canonical` database collection for the gene.
3. Otherwise prefer `CANONICAL=YES`.
4. Otherwise prefer the first `BIOTYPE=protein_coding`.
5. Otherwise use the first transcript.

The selected transcript source is stored in `INFO/selected_CSQ_criteria` as `db`, `vep`, or `random`.

### SNV data that should be avoided

Avoid these shapes because the current loader will either fail or produce poor data:

- VCFs without a `CSQ` header and `INFO/CSQ` values.
- Records with empty `CSQ` annotations.
- Missing `variant_callers`.
- Genotypes without `VAF`, `VD`, `DP`, or `GT`.
- Normal/control sample before tumor/case sample.
- Multi-allelic records if allele-specific annotations or frequencies matter. The importer joins all ALT alleles into one comma-separated string, while some frequency parsing assumes one allele string.
- Experimental transcript-only variants where all `Feature` values start with `X`, unless the CSQ gene symbol is one of `HNF1A`, `MZT2A`, `SNX9`, `KLHDC4`, `LMTK3`, or `PTPA`. Other such variants are skipped.

## CNV JSON

CNVs are loaded from a JSON object. The top-level keys are CNV identifiers, typically formatted as `chr:start-end`. The top-level key is not stored; each value object is inserted into the `cnvs` collection after the importer adds `SAMPLE_ID`.

Recommended object shape:

```json
{
  "7:55019008-55166009": {
    "nprobes": 16,
    "callers": "gatk-cnvkit",
    "ratio": 2.76289,
    "size": 147002,
    "chr": "7",
    "start": 55019008,
    "end": 55166009,
    "genes": [
      {
        "gene": "EGFR",
        "class": "somatic",
        "cnv_type": "unspec"
      }
    ]
  }
}
```

Recommended fields per CNV:

| Field | Type | Notes |
| --- | --- | --- |
| `chr` | string | Chromosome name. |
| `start` | integer | Start coordinate. |
| `end` | integer | End coordinate. |
| `size` | integer | Event size in bases. |
| `ratio` | number | Copy number/log-ratio-like value from the producing pipeline. |
| `nprobes` | integer | Number of probes/bins supporting the call. |
| `callers` | string | Caller name or combined callers, for example `gatk-cnvkit`. |
| `genes` | array | Overlapping genes. |

Recommended fields per `genes` item:

| Field | Type | Notes |
| --- | --- | --- |
| `gene` | string | Gene symbol. |
| `class` | string | Optional, for example `somatic`. |
| `cnv_type` | string | Optional, for example `gain`, `loss`, `amp`, `del`, or `unspec`. |

The loader does not require a non-empty CNV file. Empty JSON objects are accepted and result in no CNVs being inserted.

## Biomarkers JSON

Biomarkers are loaded from one JSON object. The full object is inserted into the `biomarkers` collection after the importer adds `SAMPLE_ID`.

Example:

```json
{
  "MSIS": {
    "tot": 3221,
    "som": 159,
    "perc": 4.94
  },
  "HRD": {
    "tai": 16,
    "hrd": 19,
    "lst": 6,
    "sum": 41
  },
  "name": "MadeUpSample"
}
```

Recommended conventions:

| Field | Type | Notes |
| --- | --- | --- |
| `name` | string | Sample/case name. Present in the example file and useful for traceability. |
| biomarker keys | object | Each biomarker can define its own metric object. |

The code does not enforce specific biomarker names. Current example biomarkers include `MSIS` and `HRD`.

## RNA fusions JSON

RNA fusions are loaded from a JSON array. Each item is a fusion event object. The importer adds `SAMPLE_ID` to every fusion object and inserts all objects into the `fusions` collection.

Example:

```json
[
  {
    "gene1": "RANBP2",
    "gene2": "EEF1A1",
    "genes": "RANBP2^EEF1A1",
    "calls": [
      {
        "selected": 1,
        "caller": "fusioncatcher",
        "breakpoint1": "2:108719678:+",
        "breakpoint2": "6:73520056:-",
        "spanpairs": "16",
        "spanreads": "2",
        "longestanchor": "30",
        "commonreads": "0",
        "effect": "CDS(truncated)/UTR",
        "desc": "oncogene,cancer,m11,t2,exon-exon"
      }
    ]
  }
]
```

Recommended fields per fusion:

| Field | Type | Notes |
| --- | --- | --- |
| `gene1` | string | First fusion partner. |
| `gene2` | string | Second fusion partner. |
| `genes` | string | Combined gene pair, using `^` in current examples. |
| `calls` | array | Caller-specific support entries. |

Recommended fields per `calls` item:

| Field | Type | Notes |
| --- | --- | --- |
| `selected` | integer | Optional marker for the preferred call. Current examples use `1`. |
| `caller` | string | Fusion caller name, for example `fusioncatcher`. |
| `breakpoint1` | string | Breakpoint formatted as `chrom:pos:strand`. |
| `breakpoint2` | string | Breakpoint formatted as `chrom:pos:strand`. |
| `spanpairs` | string or number | Spanning pair support. |
| `spanreads` | string or number | Spanning read support. |
| `longestanchor` | string or number | Longest anchor length. |
| `commonreads` | string or number | Common read count. |
| `effect` | string | Predicted effect, for example `in-frame` or `out-of-frame`. |
| `desc` | string | Caller/pipeline description tags. |

The current loader calls `insert_many` directly, so provide a non-empty array when loading fusions.

## Minimal YAML wiring

DNA sample with SNVs, CNVs, and biomarkers:

```yaml
groups: ["solid_GMSv3"]
name: "sample-id"
assay: "solid_GMSv3"
genome_build: 38
vcf_files: "/path/to/sample.vep.vcf.gz"
cnv: "/path/to/sample.cnv.json"
biomarkers: "/path/to/sample.biomarkers.json"
```

RNA sample with fusions:

```yaml
groups: ["rna"]
name: "sample-id-rna"
genome_build: 38
fusion_files: "/path/to/sample.fusions.json"
```

