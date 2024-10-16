from pymongo import MongoClient
from bson.objectid import ObjectId
client = MongoClient("localhost")
db = client["coyote_dev"]
var_col = db["variants_idref"]
var_col_new = db["variants"]
sample_col = db["samples"]
from collections import defaultdict
from pprint import pprint

def pick_af_fields(var):
    """
    return an AF field to top-level for variant

    Will have gnomAD_AF gnomAD genomes, then exac, then thousand genomes
    and if possible return max af for gnomad
    """
    af_dict           = {"genomad_frequency": "", "genomad_max": "", "exac_frequency": "", "1000g_frequency": ""}
    allele            = var["ALT"]
    exac              = parse_allele_freq( var["INFO"]["CSQ"][0].get("ExAC_MAF"),allele)
    thousand_g        = parse_allele_freq( var["INFO"]["CSQ"][0].get("GMAF"),allele)
    gnomad            = var["INFO"]["CSQ"][0].get("gnomAD_AF", 0)
    gnomad_genome     = var["INFO"]["CSQ"][0].get("gnomADg_AF", 0)
    gnomad_max        = var["INFO"]["CSQ"][0].get("MAX_AF", 0)

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
        af_dict["1000g_frequency"] = thousand_g
    return af_dict

def max_gnomad(gnomad):
    """
    check if gnoamd is multivalued, split and max
    """
    try:
        gnomad_list = gnomad.split('&')
        if gnomad_list:
            return float(max(gnomad_list))
    except:
        return gnomad



def parse_allele_freq(freq_str, allele):

    if freq_str:
        all_alleles = freq_str.split('&')
        for allele_frq in all_alleles:
            a = allele_frq.split(':')
            if a[0] == allele:
                return float(a[1])

    return 0

def split_on_comma(data):

    data_parts = data.split(":")
    if len(data_parts) > 1:
        return data_parts[1]
    else:
        return data

def split_on_ambersand(found_dict,string):
    """
    collect pubmed ids
    """
    try:
        string_list = string.split('&')
        for c in string_list:
            found_dict[c] = 1
        return found_dict
    except:
        found_dict[str(string)] = 1
        return found_dict

def collect_dbsnp(dbsnp_dict,dbsnp):
    """
    collect all rsids
    """
    dbsnp_list = dbsnp.split('&')
    for snp in dbsnp_list:
        if snp.startswith("rs"):
            dbsnp_dict[snp] = 1
    return dbsnp_dict
def parse_transcripts(csq:dict):
    """
    reduce redundancy from CSQ-strings
    re-organize data into top-level
    """
    transcripts    = []
    pubmed_dict    = {}
    cosmic_dict    = {}
    dbsnp_dict     = {}
    transcript_ids = {}
    hgvsc_ids      = {}
    hgvsp_ids      = {}
    gene_symbols   = {}
    for transcript in csq:
        slim_transcript = {}
        
        slim_transcript["Feature"]        = transcript.get("Feature")
        transcript_id                     = str(transcript.get("Feature")).split('.')[0]
        transcript_ids[transcript_id]     = 1
        slim_transcript["HGNC_ID"]        = transcript.get("HGNC_ID")
        gene_symbol                       = transcript.get("SYMBOL")
        slim_transcript["SYMBOL"]         = gene_symbol
        gene_symbols[gene_symbol]         = 1
        slim_transcript["PolyPhen"]       = transcript.get('PolyPhen')
        slim_transcript["SIFT"]           = transcript.get('SIFT')
        slim_transcript["Consequence"]    = transcript.get('Consequence')
        slim_transcript["ENSP"]           = transcript.get('ENSP')
        slim_transcript["BIOTYPE"]        = transcript.get('BIOTYPE')
        slim_transcript["INTRON"]         = transcript.get('INTRON')
        slim_transcript["EXON"]           = transcript.get('EXON')
        slim_transcript["CANONICAL"]      = transcript.get('CANONICAL')
        slim_transcript["MANE"]           = transcript.get('MANE_SELECT')
        slim_transcript["STRAND"]         = transcript.get('STRAND')
        slim_transcript["IMPACT"]         = transcript.get("IMPACT")

        protein_change = transcript.get('HGVSp')
        if protein_change:
            protein_change = split_on_comma(protein_change)
        slim_transcript["HGVSp"]          = protein_change
        hgvsp_ids[protein_change]         = 1
        cdna_change = transcript.get('HGVSc')
        if cdna_change:
            cdna_change = split_on_comma(cdna_change)
        slim_transcript["HGVSc"]          = cdna_change
        hgvsc_ids[cdna_change]         = 1
        cosmic = transcript.get('COSMIC')
        if cosmic:
            cosmic_dict = split_on_ambersand(cosmic_dict,cosmic)
        
        dbsnp = transcript.get("Existing_variation")
        if dbsnp:
            dbsnp_dict = collect_dbsnp(dbsnp_dict,dbsnp)

        pubmed = transcript.get("PUBMED")
        if pubmed:
            pubmed_dict = split_on_ambersand(pubmed_dict,pubmed)

        transcripts.append(slim_transcript)

    cosmic_list = list(cosmic_dict.keys())
    pubmed_list = list(pubmed_dict.keys())
    dbsnp = list(dbsnp_dict.keys())

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
    return transcripts, cosmic_list, dbsnp, pubmed_list, transcript_list_filtered, hgvsc_list_filtered, hgvsp_list_filtered, genes_list_filtered


def add_to_new_collection(var_object):
    var_col_new.insert_one(var_object)

#samples = sample_col.find( { "groups" : ["tumwgs" ]} )
#variants = (var_col.find( { "_id" : ObjectId("6706917cfcb65091269167f2") } ))
#for sample in samples:
    #sample_id = sample['_id']
variants = (var_col.find( { "SAMPLE_ID" : "67069173fcb6509126916438" } ))
for variant in variants:
    csq = variant["INFO"].get("CSQ")
    if csq:
        #print(f"variant id {variant['_id']}")
        variant.update(pick_af_fields(variant))
        slim_csq, cosmic_list, dbsnp, pubmed_list, trans, cdna, prot, genes = parse_transcripts(csq)
        variant["INFO"]["CSQ"] = slim_csq
        variant["cosmic_ids"] = cosmic_list
        variant["dbsnp_id"] = dbsnp
        variant["pubmed_ids"] = pubmed_list
        variant["HGVSp"] = prot
        variant["HGVSc"] = cdna
        variant["genes"] = genes
        variant["transcripts"] = trans
        variant["simple_id"] = f"{variant['CHROM']}_{variant['POS']}_{variant['REF']}_{variant['ALT']}"
    # add_to_new_collection(variant)
