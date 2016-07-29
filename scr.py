#filtering on phred qual scores
import variants_pb2
import sys
import os
import pysam
from pysam import VariantFile
from pymongo import MongoClient
import json
import google.protobuf.json_format as json_format
import time
import progressbar
import uuid
import google.protobuf.struct_pb2 as struct_pb2
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument("-d","--directory", help="Directory to use")
parser.add_argument("-i", "--input", help="Input file")
parser.add_argument("-db", "--database", help="Database output")
parser.add_argument("-ch", "--chromosome", help="Chromosome choice")
parser.add_argument("-pb", "--protobuf", help="Protobuf out")
p = parser.parse_args()

assert p.input

widgets = [progressbar.Timer()]
pBar = progressbar.ProgressBar(widgets=widgets, max_value=100)
vcfFile = pysam.VariantFile(p.input)
hdr = vcfFile.header
sampleNames = list(hdr.samples)
vsID = str(uuid.uuid4())
#checking if chromosome was specified
if p.chromosome is not None:
    chrom = p.chromosome
else:
    chrom = None

def main():
    if p.database is not None:
        global db
        db = get_db(p.database)
        variantset = db.VariantSet
        variantd = db.Variants
        global calls
        calls = db.Calls
        global callset
        callset = db.CallSets
    #if the directory does not already exist then create it
    if not os.path.exists("output2/variantSet/variants/calls/callsets"):
        os.makedirs("output2/variantSet/variants/calls/callsets")
    if p.protobuf is None:
        vSet_FileName = vsID + '.txt'
    else:
        vSet_FileName = vsID + '.pb'
    
    fout1 = open(os.path.join("output2/variantSet", vSet_FileName), 'w')
    if p.protobuf is None:
        fout1.write (json.dumps(json_format._MessageToJsonObject(variantSet(hdr), True)))
        fout1.close()
    else:
        fout1.write (variantSet(hdr))
        fout1.close()

    count = 0
    #checking if the database option was chosen
    if "db" in globals():
        variantset.insert_one(json_format._MessageToJsonObject(variantSet(hdr), True))
    pBar.start()
    for variant in vcfFile.fetch(chrom):
        count += 1
        if count % 100 == 0:
            pBar.update()
        if p.protobuf is None:
            v_FileName = variant.id + '.txt'
        else:
            v_FileName = variant.id + '.pb'
        #checking if the database option was chosen
        if "db" in globals():
            variantd.insert_one(json_format._MessageToJsonObject(vMes(variant), True))
        if p.protobuf is None:
            if not os.path.isfile(v_FileName):
                fout2 = open(os.path.join("output2/variantSet/variants", v_FileName), 'w')
            fout2.write (json.dumps(json_format._MessageToJsonObject(vMes(variant))))
        else:
            if not os.path.isfile(v_FileName):
                fout2 = open(os.path.join("output2/variantSet/variants", v_FileName), 'w')
            fout2.write (vMes(variant))
    pBar.finish()
    fout2.close()  

def get_db(pdatabase):
    client = MongoClient() 
    db = client[pdatabase]
    return db

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
    if p.protobuf is not None:
        return gaVariantVS.SerializeToString()
    else:
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
    if p.protobuf is None:
        cs_txt_FileName = cs_FileName + '.txt'
    else:
        cs_txt_FileName = cs_FileName + '.pb'
    if not os.path.isfile(cs_txt_FileName):
            fout4 = open(os.path.join("output2/variantSet/variants/calls/callsets", cs_txt_FileName), 'w')
    if p.protobuf is None:
        fout4.write (json.dumps(json_format._MessageToJsonObject(gaVariantCS, True)))
    else:
        fout4.write (gaVariantCS.SerializeToString())
    fout4.close()
    if "db" in globals():
        callset.insert_one(json_format._MessageToJsonObject(gaVariantCS, True))
    return gaVariantCS

def callMes(call_record, sample_name, variant_id):
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
    if variant_id is not None:
        gaVariantC.info["variant_id"].append(variant_id)
    c_FileName = gaVariantC.call_set_id
    if p.protobuf is None:
        c_txt_FileName = c_FileName + '.txt'
    else:
        c_txt_FileName = c_FileName + '.pb'
    if not os.path.isfile(c_txt_FileName):
            fout3 = open(os.path.join("output2/variantSet/variants/calls", c_txt_FileName), 'w')
    if p.protobuf is None:    
        fout3.write (json.dumps(json_format._MessageToJsonObject(gaVariantC, True)))
    else:
        fout3.write (gaVariantC.SerializeToString())
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
    if variant.id is not None:
        gaVariant.names.append(variant.id)
        for sample_name in sampleNames:
            call_record = variant.samples[sample_name]
            gaVariant.calls.extend([callMes(call_record,sample_name,variant.id)])
    if p.protobuf is not None:
        return gaVariant.SerializeToString()
    else:
        return gaVariant

main()