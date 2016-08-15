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
        global variantset
        variantset = db.VariantSet
        global variantd
        variantd = db.Variants
        global calls
        calls = db.Calls
        global callset
        callset = db.CallSets
#if the directory does not already exist then create it
    if not os.path.exists("output2/variantSet/variants/calls/callsets"):
        os.makedirs("output2/variantSet/variants/calls/callsets")
    variantSet(hdr)
#count used to update the progressbar
    count = 0
    pBar.start()
    for variant in vcfFile.fetch(chrom):
        count += 1
        vMes(variant)
        if count % 100 == 0:
            pBar.update()
    pBar.finish()        
def get_db(pdatabase):
    client = MongoClient() 
    db = client[pdatabase]
    return db

"""def test_connection(client):
    if client.alive() //alive func no longer exists == False:
        return "no mongodb connection"
    else:
        return"""
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
#Contains VCF Header information and calls the vsMetadata function.
def vHeader(hdr):
    formats = hdr.formats.items()
    infos = hdr.info.items()
    meta = []
    for prefix, content in [("Format", formats), ("Info", infos)]:
        for key, value in content:
            meta.append(vsMetadata(key,value.type,value.number,value.description))
    return meta     
#Creates GA4GH Variant Set message.
#Also contains the files metadata by calling vHeader function.
def variantSet(hdr):
    ranId = uuid.uuid4()
    gaVariantVS = variants_pb2.VariantSet()
    gaVariantVS.reference_set_id = hdr.rid
    gaVariantVS.id = str(vsID)
    gaVariantVS.name = str(hdr.contigs)
    gaVariantVS.dataset_id = str(ranId)
    gaVariantVS.metadata.extend(vHeader(hdr))
    if p.protobuf is None:
        vSet_FileName = vsID + '.txt'
    else:
        vSet_FileName = vsID + '.pb'
    fout1 = open(os.path.join("output2/variantSet", vSet_FileName), 'w')
    if p.protobuf is None:
        fout1.write (json.dumps(json_format._MessageToJsonObject(gaVariantVS, True)))
    else:
        fout1.write (gaVariantVS.SerializeToString())
    fout1.close()
    if "db" in globals():
        variantset.insert_one(json_format._MessageToJsonObject(gaVariantVS, True))
#Creates a GA4GH Call Set Message
def callSet(sampleNames, callSetId):
    gaVariantCS = variants_pb2.CallSet()
    gaVariantCS.name = str(sampleNames)
    #gaVariantCS.bio_sample_id = //Leave blank
    gaVariantCS.variant_set_ids.append(str(vsID))
    gaVariantCS.created = int(time.time())
    gaVariantCS.updated = int(time.time())
    #gaVariantCS.info = //Not currently utilized
    cs_FileName = callSetId
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
    return
#Creates a GA4GH Call message.
#Calls the CallSet function, sends it sampleNames and the call set ID.
def callMes(call_record, sample_name, variant_id):
    gaVariantC = variants_pb2.Call()
    call_set_id = uuid.uuid4()
    gaVariantC.call_set_name = sample_name
    gaVariantC.call_set_id = str(call_set_id)
    try:
        if call_record.allele_indices is not None:
            gaVariantC.genotype.extend(list(call_record.allele_indices))
    except:
        pass
    if call_record.phased:
        phaseset = str(call_record.phased)  
    gaVariantC.phaseset = str(phaseset)
    for key, value in call_record.iteritems():
        if key == 'GL' and value is not None:
            gtlikelihood = value
            gaVariantC.genotype_likelihood.extend(list(gtlikelihood)) #GTLikelihood is not always in a VCF File
    gaVariantC.info["variant_id"].append(str(variant_id))
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
    callSet(sampleNames, gaVariantC.call_set_id)
    return gaVariantC

#Creates a GA4GH variant message.  
#Calls the Call Message function also.
def vMes(variant):
    ranId = uuid.uuid4()
    gaVariant = variants_pb2.Variant()
    gaVariant.variant_set_id = str(vsID)
    variant_id = ranId
    gaVariant.id = str(variant_id)
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
    else:
        variant.id = None
    for sample_name in sampleNames:
        call_record = variant.samples[sample_name]
        gaVariant.calls.extend([callMes(call_record,sample_name,variant_id)])
    if p.protobuf is None:
        v_FileName = gaVariant.id + '.txt'
    else:
        v_FileName = gaVariant.id + '.pb'
    if not os.path.isfile(v_FileName):
        fout2 = open(os.path.join("output2/variantSet/variants", v_FileName), 'w')
    if "db" in globals():
        variantd.insert_one(json_format._MessageToJsonObject(gaVariant, True))
    if p.protobuf is None:
        fout2.write (json.dumps(json_format._MessageToJsonObject(gaVariant, True)))
        fout2.close()
    else:
        fout2.write (gaVariant.SerializeToString())
        fout2.close() 
    return

main()
