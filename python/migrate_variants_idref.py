from pymongo import MongoClient
from bson.objectid import ObjectId

client = MongoClient("localhost")
db = client["coyote_dev_3"]
var_col = db["variants_idref"]
var_col_new = db["variants"]
var_col_filtered = db["variants_filtered"]
sample_col = db["samples"]
canonical_col = db["refseq_canonical"]

from collections import defaultdict
from pprint import pprint

##### TODO: TO SELECT CSQ (MODIFY FOT THE NEW DATA FROM NEW VEP) #####
## FUNCTIONS RELATED FOR SELECTING CSQ ##


def get_canonical():
    canonical_dict = {}
    canonical = canonical_col.find({})
    for c in canonical:
        canonical_dict[c["gene"]] = c["canonical"]

    return canonical_dict


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


##### END OF FUNCTIONS RELATED FOR SELECTING CSQ #####


def pick_af_fields(var):
    """
    return an AF field to top-level for variant

    Will have gnomAD_AF gnomAD genomes, then exac, then thousand genomes
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

        # removes if the transcript is experimental and not in the list of genes
        if transcript.get("Feature").startswith("X") and transcript.get(
            "SYMBOL"
        ) not in ["HNF1A", "MZT2A", "SNX9", "KLHDC4", "LMTK3", "PTPA"]:
            continue

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


def add_to_new_collection(var_object):
    var_col_new.insert_one(var_object)


def add_to_filtered_collection(var_object):
    var_col_filtered.insert_one(var_object)


## CANONICAL TRANSCRIPTS
canonical_dict = get_canonical()

# samples = sample_col.find( { "groups" : ["tumwgs" ]} )
# variants = (var_col.find( { "_id" : ObjectId("6706917cfcb65091269167f2") } ))
# for sample in samples:
# sample_id = sample['_id']
# variants = (var_col.find( { "SAMPLE_ID" : "67069173fcb6509126916438" } ))
# variants = (var_col.find( { "SAMPLE_ID" : "6618fd85dbdd45ba0b732513" } ))
# variants = (var_col.find( { "SAMPLE_ID" : "65fbdec69bc47223004e97e1" } )) # solid_integration_test-2
# variants = (var_col.find( { "SAMPLE_ID" : "659b37c09bc4724425065161" } )) # 23MD14600-wgs
# variants = (var_col.find( { "SAMPLE_ID" : "66a0a4f09bc472c8ed4acfe1" } )) # 23PL18971-00-01-02-DNA

for sample_id in [
    "6706917cfcb65091269167f2",
    "67069173fcb6509126916438",
    "6618fd85dbdd45ba0b732513",
    "65fbdec69bc47223004e97e1",
    "659b37c09bc4724425065161",
    "66a0a4f09bc472c8ed4acfe1",
]:
    variants = var_col.find({"SAMPLE_ID": sample_id})

    for variant in variants:
        csq = variant["INFO"].get("CSQ")
        # removes the variant all together if the csq only has experimental transcripts but retains if the any of those experimental transcripts are in the list of genes
        valid_var = True
        if csq:
            all_features = [c.get("Feature") for c in csq]
            all_X_genes = [
                c.get("SYMBOL") for c in csq if c.get("Feature", "").startswith("X")
            ]

            if all([f.startswith("X") for f in all_features]) and not any(
                [
                    g in ["HNF1A", "MZT2A", "SNX9", "KLHDC4", "LMTK3", "PTPA"]
                    for g in list(set(all_X_genes))
                ]
            ):
                valid_var = False

        if valid_var:
            # print(f"variant id {variant['_id']}")
            variant.update(pick_af_fields(variant))
            (
                slim_csq,
                cosmic_list,
                dbsnp,
                pubmed_list,
                trans,
                cdna,
                prot,
                genes,
                hotspots_dict,
            ) = parse_transcripts(csq)
            selected_csq, selected_csq_source = select_csq(
                slim_csq, canonical_dict
            )  # How slow is this

            variant["INFO"]["CSQ"] = selcted_transcript_removal(
                slim_csq, selected_csq["Feature"]
            )
            variant["INFO"]["selected_CSQ"] = selected_csq
            variant["INFO"]["selected_CSQ_criteria"] = selected_csq_source
            variant["selected_csq_feature"] = selected_csq["Feature"]
            variant["variant_class"] = selected_csq.get("VARIANT_CLASS")
            variant["cosmic_ids"] = cosmic_list
            variant["dbsnp_id"] = dbsnp
            variant["pubmed_ids"] = pubmed_list
            variant["HGVSp"] = prot
            variant["HGVSc"] = cdna
            variant["genes"] = genes
            variant["transcripts"] = trans
            variant["hotspots"] = [hotspots_dict]
            variant["simple_id"] = (
                f"{variant['CHROM']}_{variant['POS']}_{variant['REF']}_{variant['ALT']}"
            )
            # add_to_new_collection(variant)
        else:
            # Adds all the removed variants to a new collection for further inspection
            # add_to_filtered_collection(variant)
            pass
