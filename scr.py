import variants_pb2
import sys
import os
import pysam
from pysam import VariantFile

import json
import google.protobuf.json_format as json_format
import time
import progressbar
import uuid
import google.protobuf.struct_pb2 as struct_pb2

#for command line arguments
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d","--directory", help="directory to use")
parser.add_argument("-i", "--input", help="Input file")
parser.add_argument("-o", "--output", help="Output file")
p = parser.parse_args()

assert p.directory
assert p.input and p.output
#progress indicator for file
widgets = [progressbar.Timer()]
pBar = progressbar.ProgressBar(widgets=widgets, max_value=100).start()

#Main
vcfFile = pysam.VariantFile(p.input)
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


def vsMetadata(key, type1, number, description):
	gaVariant_metaData = variants_pb2.VariantSetMetadata()
	gaVariant_metaData.key = key
	#gaVariant_metaData.value = value //same as key
	gaVariant_metaData.type = type1
	gaVariant_metaData.number = str(number)
	gaVariant_metaData.description = description
	return gaVariant_metaData

def vHeader(hdr):
 	formats = hdr.formats.items()
	infos = hdr.info.items()
	meta = []
	for prefix, content in [("Format", formats), ("Info", infos)]:
		for key, value in content:
			meta.append(vsMetadata(key,value.type,value.number,value.description))
	return meta		

def variantSet(hdr):
	ranId = uuid.uuid4()
	gaVariantVS = variants_pb2.VariantSet()
	#gaVariantVS.reference_set_id = 
	gaVariantVS.id = str(vsID)
	gaVariantVS.name = list(str(hdr.contigs))
	gaVariantVS.dataset_id = str(ranId)
	gaVariantVS.metadata.extend(vHeader(hdr))
	return gaVariantVS

def callSet():
	gaVariantCS = variants_pb2.CallSet()
	gaVariantCS.name = str(sampleNames)
	#gaVariantCS.bio_sample_id = //Leave blank
	gaVariantCS.variant_set_ids.append(str(vsID))
	gaVariantCS.created = int(time.time())
	gaVariantCS.updated = int(time.time())
	#gaVariantCS.info = //seems useless
	return gaVariantCS

def callMes(call_record, sample_name):
	gaVariantC = variants_pb2.Call()
	call_set_id = uuid.uuid4()
	gaVariantC.call_set_name = sample_name
	gaVariantC.call_set_id = str(call_set_id)
	gaVariantC.genotype.extend(list(call_record.allele_indices))
	if call_record.phased:
		phaseset = str(call_record.phased)	
	gaVariantC.phaseset = str(phaseset)
	for key, value in call_record.iteritems():
		if key == 'GL' and value is not None:
			gtlikelihood = value
			gaVariantC.genotype_likelihood.extend(list(gtlikelihood))
	#gaVariantC.info 
	callSet()
	return gaVariantC

def vMes(variant):
	ranId = uuid.uuid4()
	gaVariant = variants_pb2.Variant()
	gaVariant.variant_set_id = str(vsID)
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
	for sample_name in sampleNames:
		call_record = variant.samples[sample_name]
		gaVariant.calls.extend([callMes(call_record,sample_name)])
	return gaVariant

fout = open(p.output,"w")
fout.write (json.dumps(json_format._MessageToJsonObject(variantSet(hdr), True)))
for variant in vcfFile.fetch(chromo, 0, 100):
	fout.write (json.dumps(json_format._MessageToJsonObject(vMes(variant), True)))
fout.close()