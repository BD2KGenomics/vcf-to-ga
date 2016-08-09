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
var_ann_set_id = uuid.uuid4()
def main():


def AnalysisRes():
	analysis_id
	result
	score
#def AnalysisMes():

def HGVSann():
	genomic
	transcript
	protein
def alleleLoc():
	start
	end
	reference_sequence
	alternate_sequence
def TranscEff():
	id
	feature_id
	alternate_bases
	(repeatd OntologyTerm) effects
	hgvs_annotation
	cdna_location
	protein_location
	(repeated AnalysisResult) analysis_result

def VarAnnSet(gaVariantVS_id):
	vaSet = allele_annotations_pb2.VariantAnnotationSet()
	vaSet.id = str(var_ann_set_id)
	vaSet.variant_set_id = gaVariantVS_id
	vaSet.name =
	vaSet.analysis =
def VarAnnMes(variant_record, gaVariant_id):
	vAnMes = allele_annotations_pb2.VariantAnnotation()
	ranId = uuid.uuid4()
	vAnMes.id = str(ranId)
	vAnMes.variant_id = gaVariant_id
	vAnMes.variant_annotation_set_id
	vAnMes.created
	(repeated Transcript Effect) vAnMes.transcript_effects
	vAnMes.info

