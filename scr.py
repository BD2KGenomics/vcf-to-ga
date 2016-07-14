import variants_pb2
import sys
import os
import pysam
from pysam import VariantFile

import google.protobuf.json_format as json_format
import time
import uuid
import google.protobuf.struct_pb2 as struct_pb2

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", help="Input file")
parser.add_argument("-o", "--output", dest="output", help="Output file")
p = parser.parse_args()

assert p.input and p.output


#Main
vcfFile = pysam.VariantFile(p.input)
sys.stdout = open(p.output,"w")
hdr = vcfFile.header
#chromosome choice
chromo = "ref_brca1"
sampleNames = list(hdr.samples)
vsID = uuid.uuid4()
#this function taken from ga4gh/datamodel/variants.py.
def _encodeValue(value):
    if isinstance(value, (list, tuple)):
        return [struct_pb2.Value(string_value=str(v)) for v in value]
    else:
        return [struct_pb2.Value(string_value=str(value))]


def vHeader(hdr):
	ranId = uuid.uuid4()
 	gaVariantMD = variants_pb2.VariantSetMetadata()
 	gaVariantMD.info = hdr.info.items()
 	gaVariantMD.id = str(ranId)
 	formats = hdr.formats.items()
 	gaVariantMD.type = str(formats).strip('[]')
 	#gaVariantMD.number = 
 	gaVariantMD.description = str(vcfFile.description)
 	#gaVariantMD.key = str(list(hdr.info))
 	#gaVariantMD.value = 
 	return gaVariantMD

def variantSet(variant):
	ranId = uuid.uuid4()
	gaVariantVS = variants_pb2.VariantSet()
	gaVariantVS.reference_set_id = str(variant.rid)
	gaVariantVS.id = str(vsID)
	gaVariantVS.name = variant.contig
	gaVariantVS.dataset_id = str(ranId)
	#gaVariantVS.metadata = 
	return gaVariantVS

def callSet(variant):
	gaVariantCS = variants_pb2.CallSet()
	ranId = uuid.uuid4()
	gaVariantCS.name = str(sampleNames)
	gaVariantCS.bio_sample_id = str(ranId)
	gaVariantCS.variant_set_ids.append(str(vsID))
	gaVariantCS.created = int(time.time())
	gaVariantCS.updated = int(time.time())
	#gaVariantCS.info map
	return gaVariantCS

def callMes(variant):
	ranId = uuid.uuid4()
	gaVariantC = variants_pb2.Call()
	cs = callSet(variant)
	gaVariantC.call_set_name = cs.name
	gaVariantC.call_set_id = str(ranId)
	#gaVariantC.genotype = where are you
	#gaVariantC.phaseset =
	#for key, value in 
	#gaVariant.genotype_likelihood
	#gaVariantC.info = #optional
	return gaVariantC

def vMes(variant):
	ranId = uuid.uuid4()
	gaVariant = variants_pb2.Variant()
	call = variantSet(variant)
	gaVariant.variant_set_id = call.id
	gaVariant.id = str(ranId)
	gaVariant.reference_name = variant.contig
	if variant.id is not None:
		gaVariant.names.append(variant.id)
	gaVariant.created = int(time.time())
	gaVariant.updated = int(time.time())
	gaVariant.start = variant.start
	gaVariant.end = variant.stop
	gaVariant.reference_bases = variant.ref
	if variant.alts is not None:
		gaVariant.alternate_bases.extend(list(variant.alts))
	for key, value in variant.info.iteritems():
		if value is not None:
			gaVariant.info[key].values.extend(_encodeValue(value))
	#for calls in callMes(variant):
		#gaVariant.calls.add().
	return gaVariant


#print (json_format._MessageToJsonObject(vHeader(hdr), True))
for variant in vcfFile.fetch(chromo, 0, 100):
	print (json_format._MessageToJsonObject(vMes(variant), True))