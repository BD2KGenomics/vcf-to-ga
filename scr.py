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
def vHeader(hdr):
	ranId = uuid.uuid4()
 	gaVariantMD = variants_pb2.VariantSetMetadata()
 	for key,value in hdr.info.iteritems():
		gaVariantMD.info[key]
 	gaVariantMD.id = str(ranId)
 	#gaVariantMD.type = hdr.
 	#gaVariantMD.number = 
 	gaVariantMD.description = str(vcfFile.description)
 	gaVariantMD.key = str(_encodeValue(hdr.info.keys))
 	for value in hdr.info.itervalues():
 		#if value is not None:
 		gaVariantMD.value = str(value)
 	return gaVariantMD

def vsMes(rec):
	ranId = uuid.uuid4()
	gaVariantVS = variants_pb2.VariantSet()
	gaVariantVS.reference_set_id = str(ranId) #pysam has an option of rec.rid (reference id, but shows up as just a 0)
	gaVariantVS.id = str(ranId)
	gaVariantVS.name = rec.contig
	gaVariantVS.dataset_id = str(ranId)
	#gaVariantVS.metadata = 
	return gaVariantVS

def csMes(hdr, rec):
	ranId = uuid.uuid4()
	gaVariantCS = variants_pb2.CallSet()
	gaVariantCS.name = str(list(hdr.samples))
	gaVariantCS.bio_sample_id = str(ranId)
	#gaVariantCS.variant_set_ids = gaVariantVS.id
	gaVariantCS.created = int(time.time())
	gaVariantCS.updated = int(time.time())
	#gaVariantCS.info map
	return gaVariantCS

def callMes(rec,hdr):
	ranId = uuid.uuid4()
	gaVariantC = variants_pb2.Call()
	#gaVariantC.call_set_name = gaVariantCS.name
	gaVariantC.call_set_id = str(ranId)
	#for gt,rec in v.samples.items():
		#gaVariantC.genotype =
	#gaVariantC.phaseset =
	#for key,value in rec.info:
		#if value is not None:
			#gaVariantC.genotype_likelihood
	#gaVariantC.info = #optional
	return gaVariantC

def vMes(rec,hdr):
	ranId = uuid.uuid4()
	gaVariant = variants_pb2.Variant()
	#gaVariant.variant_set_id =
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
			gaVariant.info[key].values.extend(_encodeValue(value))
	#for callsetid in callsetids:
		#gaVariant.calls.add()
	return gaVariant

#Main 
#output raw protobuf
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
for rec in vcfFile.fetch(chromo, 0, 100):
	print (json_format._MessageToJsonObject(vMes(rec,hdr), True))
	print (json_format._MessageToJsonObject(vsMes(rec), True))
	print (json_format._MessageToJsonObject(csMes(hdr,rec), True))
	print (json_format._MessageToJsonObject(callMes(rec,hdr), True))
