import variants_pb2
import sys
import pysam
from pysam import VariantFile
import google.protobuf.json_format as json_format

def toga(rec):
	gaVariant = variants_pb2.variants()
	gavariant.variant = rec
	


#VCF file to read
file = sys.argv[1]
#file to output to
ofile = sys.argv[2]
sys.stdout = open(ofile, "w")
vcf = pysam.VariantFile(file)
#chromosome choice
chromo = "ref_brca1"
for rec in vcf.fetch(chromo, 0, 100):
	print (json_format._MessageToJsonObject(varmes(rec), True))