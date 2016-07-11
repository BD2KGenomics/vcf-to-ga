import variants_pb2
import sys

import pysam
from pysam import VariantFile

import google.protobuf.json_format as json_format
import time
import uuid
import google.protobuf.struct_pb2 as struct_pb2

#this function taken from ga4gh/datamodel/variants.py.
def _encodeValue(value):
    if isinstance(value, (list, tuple)):
        return [struct_pb2.Value(string_value=str(v)) for v in value]
    else:
        return [struct_pb2.Value(string_value=str(value))]

def vMes(rec,hdr):
	ranId = uuid.uuid4()
	gaVariant = variants_pb2.Variant()
	gaVariant.id = str(ranId)
	gaVariant.reference_name = rec.contig
	if rec.id is not None:
		gaVariant.names.append(rec.id)
	gaVariant.created = int(time.time())
	gaVariant.updated = int(time.time())
	gaVariant.start = rec.start
	gaVariant.end = rec.stop
	gaVariant.reference_bases = rec.ref
	if rec.alts is not None:
		gaVariant.alternate_bases.extend(list(rec.alts))
	for key, value in rec.info.iteritems():
		if value is not None:
			gaVariant.info[key].values.extend(_encodeValue(value)) #gives typeError without encodeValue function
	return gaVariant

def vHeader(hdr):
 	gaVariantMD = variants_pb2.VariantSetMetadata()
 	#for key, value in hdr.info.iteritems():
 		#if value is not None:
 			#gaVariantMD.info[key].values.extend()

 	#gaVariantMD.value = hdr.info.values
 	#gaVariantMD.id = 
 	#gaVariantMD.type = 
 	#gaVariantMD.number = hdr.value
 	gaVariantMD.description = str(vcfFile.description)
 	return gaVariantMD

def vsMes(rec):
	gaVariantVS = variants_pb2.VariantSet()
	gaVariantVS.reference_set_id = str(rec.rid)
	#gaVariantVS.id = 
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
for rec in vcfFile.fetch(chromo, 0, 100):
	print (json_format._MessageToJsonObject(vMes(rec,hdr), True))
	print (json_format._MessageToJsonObject(vsMes(rec), True))