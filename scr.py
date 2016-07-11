import variants_pb2
import sys
import pysam
from pysam import VariantFile
import google.protobuf.json_format as json_format

def vMes(rec,hdr):
	gaVariant = variants_pb2.Variant()
	gaVariant.id = rec.id
	gaVariant.reference_name = rec.chrom
	#gaVariant.variant_set_id = 
	#gaVariant.names =
	#gaVariant.created = 
	#gaVariant.updated = 
	gaVariant.start = rec.start
	gaVariant.end = rec.stop
	gaVariant.reference_bases = rec.ref
	if rec.alts is not None:
		gaVariant.alternate_bases.extend(list(rec.alts))

	"""gaVariant. = hdr.info
 	gaVariant. = hdr.version
	gaVariant. = hdr.samples
 	gaVariant. = hdr.records
 	gaVariant. = hdr.contigs
 	gaVariant. = hdr.filters
 	gaVariant. = hdr.formats"""
 	return gaVariant

def vHeader(hdr):
 	gaVariantMD = variants_pb2.VariantSetMetadata()
 	#for keys in hdr.info.keys:
 		#gaVariantMD.key
 	#gaVariantMD.value = hdr.info.values
 	gaVariantMD.id = str(hdr.info.items)
 	#gaVariantMD.type = 
 	#gaVariantMD.number = hdr.value
 	gaVariantMD.description = vcfFile.description
 	return gaVariantMD

def vsMes(rec):
	gaVariantVS = variants_pb2.VariantSet()
	gaVariantVS.reference_set_id = str(rec.rid)
	gaVariantVS.id = rec.id
	#gaVariantVS.name = 
	#gaVariantVS.dataset_id = 
	#gaVariantVS.metadata =
	return gaVariantVS
#Main 

#VCF file to read
file = sys.argv[1]
#file to output to
ofile = sys.argv[2]
sys.stdout = open(ofile, "w")
vcfFile = pysam.VariantFile(file)
hdr = vcfFile.header
#chromosome choice
chromo = "ref_brca1"
print (json_format._MessageToJsonObject(vHeader(hdr), True))
for rec in vcfFile.fetch(chromo, 0, 1000):
	print (json_format._MessageToJsonObject(vMes(rec,hdr), True))
	print (json_format._MessageToJsonObject(vsMes(rec), True))
	
