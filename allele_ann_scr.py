import allele_annotations_pb2
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

def main():
	for variant_record in vcfFile.fetch():
		VarAnnMes(variant_record)

def AnalysisRes():

def AnalysisMes():

def TranscEff():

def VarAnnSet():

def VarAnnMes(variant_record):
	vAnMes = allele_annotations_pb2.VariantAnnotation()
	for annotations in variant_record:

	vAnMes.id
	vAnMes.variant_id
	vAnMes.variant_annotation_set_id
	vAnMes.created
	(repeated Transcript Effect) vAnMes.transcript_effects
	vAnMes.info

