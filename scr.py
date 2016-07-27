import variants_pb2
import sys
import os
import pysam
from pysam import VariantFile
from pymongo import MongoClient
import json
import google.protobuf.json_format as json_format
import time
import datetime
import progressbar
import uuid
import google.protobuf.struct_pb2 as struct_pb2

#command line argument parser
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d","--directory", help="directory to use")
parser.add_argument("-i", "--input", help="Input file")
parser.add_argument("-o", "--output", help="Output file")
parser.add_argument("-db", "--database", help="database output")
parser.add_argument("-ch", "--chromosome", help="Chromosome choice")
p = parser.parse_args()

assert p.input
#assert p.chromosome

#progress indicator for file need to add update
#widgets = [progressbar.Timer()]
#pBar = progressbar.ProgressBar(widgets=widgets, max_value=100)

def get_db(pdatabase):
	client = MongoClient()
	db = client.test_vcf
	return db

vcfFile = pysam.VariantFile(p.input)
hdr = vcfFile.header
sampleNames = list(hdr.samples)
vsID = str(uuid.uuid4())

if p.chromosome is not None:
	chrom = p.chromosome
else:
	chrom = "ref_brca1"

if p.database is not None:
	db = get_db(p.database)
	variantset = db.VariantSet
	variantd = db.Variants
	calls = db.Calls
	callset = db.CallSets
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
	gaVariantVS.name = str(hdr.contigs)
	gaVariantVS.dataset_id = str(ranId)
	gaVariantVS.metadata.extend(vHeader(hdr))
	return gaVariantVS

def callSet(sampleNames):
	gaVariantCS = variants_pb2.CallSet()
	gaVariantCS.name = str(sampleNames)
	#gaVariantCS.bio_sample_id = //Leave blank
	gaVariantCS.variant_set_ids.append(str(vsID))
	gaVariantCS.created = int(time.time())
	gaVariantCS.updated = int(time.time())
	#gaVariantCS.info = //seems useless
	cs_FileName = "Call_Sets"
	cs_txt_FileName = cs_FileName + '.txt'
	if not os.path.isfile(cs_txt_FileName):
			fout4 = open(os.path.join("output2/variantSet/variants/calls/callsets", cs_txt_FileName), 'w')
			fout4.write (json.dumps(json_format._MessageToJsonObject(gaVariantCS, True)))
	fout4.close()
	if "db" in globals():
		callset.insert_one(json_format._MessageToJsonObject(gaVariantCS, True))
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
	c_FileName = gaVariantC.call_set_id
	c_txt_FileName = c_FileName + '.txt'
	if not os.path.isfile(c_txt_FileName):
			fout3 = open(os.path.join("output2/variantSet/variants/calls", c_txt_FileName), 'w')
			fout3.write (json.dumps(json_format._MessageToJsonObject(gaVariantC, True)))
	fout3.close()
	if "db" in globals():
		calls.insert_one(json_format._MessageToJsonObject(gaVariantC, True))
	callSet(sampleNames)
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

if not os.path.exists("output2/variantSet/variants/calls/callsets"):
	os.makedirs("output2/variantSet/variants/calls/callsets")
vSet_FileName = vsID + '.txt'
fout1 = open(os.path.join("output2/variantSet", vSet_FileName), 'w')
fout1.write (json.dumps(json_format._MessageToJsonObject(variantSet(hdr), True)))
fout1.close()

if "db" in globals():
	variantset.insert_one(json_format._MessageToJsonObject(variantSet(hdr), True))

#pBar.start()
for variant in vcfFile.fetch(chrom):
	v_FileName = variant.id + '.txt'
	if "db" in globals():
		variantd.insert_one(json_format._MessageToJsonObject(vMes(variant), True))
	if not os.path.isfile(v_FileName):
		fout2 = open(os.path.join("output2/variantSet/variants", v_FileName), 'w')
		time_begin = datetime.datetime.fromtimestamp(time.time())
		fout2.write (json.dumps(json_format._MessageToJsonObject(vMes(variant), True)))
time_end = datetime.datetime.fromtimestamp(time.time())		
fout2.close()
print("Time elapsed: ", str(time_end - time_begin))