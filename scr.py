import variants_pb2
import sys
import pysam
from pysam import VariantFile
import google.protobuf.json_format as json_format
import datetime


def vMes(rec,hdr):
	gaVariant = variants_pb2.Variant()
	now = datetime.datetime.now()
	gaVariant.id = rec.id
	gaVariant.reference_name = rec.contig
	#gaVariant.names =
	gaVariant.created = int(now.microsecond)
	gaVariant.updated = int(now.microsecond)
	gaVariant.start = rec.start
	gaVariant.end = rec.stop
	gaVariant.reference_bases = rec.ref
	if rec.alts is not None:
		gaVariant.alternate_bases.extend(list(rec.alts))
	for key, value in rec.info.iteritems():
		if value is not None:
			gaVariant.info[key]

 	return gaVariant

def vHeader(hdr):
 	gaVariantMD = variants_pb2.VariantSetMetadata()
 	for key, value in hdr.info.iteritems():
 		if value is not None:
 			gaVariantMD.info[key]

 	#gaVariantMD.value = hdr.info.values
 	#gaVariantMD.id = 
 	#gaVariantMD.type = 
 	#gaVariantMD.number = hdr.value
 	gaVariantMD.description = str(vcfFile.description)
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
#output raw protobuf
#VCF file to read
#make note of fields not in file
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