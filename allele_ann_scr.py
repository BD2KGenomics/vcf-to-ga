import allele_annotations_pb2
from scr import variantSet, vMes
import sys
import os
import pysam
from pysam import VariantFile
from pymongo import MongoClient
import json
import google.protobuf.json_format as json_format
import time
import uuid
import google.protobuf.struct_pb2 as struct_pb2
import argparse

#This function creates analysis result object
def AnalysisRes():
    aRes = allele_annotations_pb2.AnalysisResult()
    aRes.analysis_id
    aRes.result
    aRes.score

#this function creates allele location object
def alleleLoc():
    aLoc = allele_annotations_pb2.AlleleLocation()
    aLoc.start
    aLoc.end
    aLoc.reference_sequence
    aLoc.alternate_sequence

#this function creates variant annotation set object
def VarAnnSet(gaVariantVS_id):
    vAnSet = allele_annotations_pb2.VariantAnnotationSet()
    vAnSet.id = var_ann_set_id
    vAnSet.variant_set_id = gaVariantVS_id
    vAnSet.name =
    vAnSet.analysis =

#This function creates hgvs annotation object
def hgvsAnn():
    hgvs = allele_annotations_pb2.HGVSAnnotation()
    hgvs.genomic
    hgvs.transcript
    hgvs.protein
    return hgvs

#This function creates transcript effect object
def TranscEff(ann):
    tEff = allele_annotations_pb2.TranscriptEffect()
    tEff.id =
    tEff.feature_id =
    tEff.alternate_bases =
    (repeatd OntologyTerm) tEff.effects =
    tEff.hgvs_annotation.extend(hgvsAnn()) =
    tEff.cdna_location =
    tEff.protein_location =
    tEff.analysis_result.extend(AnalysisRes()) =
    return tEff

#this function creates variant annotation message object
def VarAnnMes(variant_record, gaVariant_id):
    vAnMes = allele_annotations_pb2.VariantAnnotation()
    ranId = uuid.uuid4()
    vAnMes.id = str(ranId)
    vAnMes.variant_id = gaVariant_id
    vAnMes.variant_annotation_set_id = var_ann_set_id
    vAnMes.created = int(time.time())
    for ann in rec.info["ANN"]: #going to turn into function and send Type[x] 
        Type = ann.split("|")
        allele = Type[0]
        annotation = Type[1]
        ann_impact = Type[2]
        Gene_Name = Type[3]
        Gene_ID = Type[4]
        Feature_Type = Type[5]
        Feature_ID = Type[6]
        Transcript_BioType = Type[7]
        Rank = Type[8]
        HGVSc = Type[9]
        HGVSp = Type[10]
        cDNA = Type[11]
        CDS = Type[12]
        AA = Type[13]
        Distance = Type[14]
        EWI = Type[15]
    vAnMes.info

#sets command line input file 
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file")
p = parser.parse_args()
assert p.input

vcfFile = pysam.VariantFile(p.input)
hdr = vcfFile.header
var_ann_set_id = str(uuid.uuid4())

for variant_record in vcfFile.fetch():
    VarAnnMes(variant_record, gaVariant.id)