from pysam import VariantFile
from pprint import pprint
import cmdvcf
import sys

infile = sys.argv[1]
vcf_object = VariantFile(infile)

count = 0
# stream one variant at the time, do what you wish with it. Either CMD-style or pysam style
for var in vcf_object.fetch():
    var_dict = cmdvcf.parse_variant(var,vcf_object.header)
    #pprint(var_dict)
    if count == 5:
        exit()
    count+=1

## loads full vcf into memory, variants is a list of var-dicts
#header,variants = cmdvcf.parse_vcf(bcf_in)

