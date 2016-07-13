import variants_pb2
import sys
import pysam
from pysam import VariantFile
import google.protobuf.json_format as json_format

def vMes(rec):
	gaVariant = variants_pb2.Variant()
	gaVariant.id = rec.id
	gaVariant.reference_name = rec.chrom
	gaVariant.variant_set_id = 
	#gaVariant.names =
	#gaVariant.created = 
	#gaVariant.updated = 
	gaVariant.start = rec.start
	gaVariant.end = rec.stop
	gaVariant.reference_bases = rec.ref
	if rec.alts is not None:
		gaVariant.alternate_bases.extend(list(rec.alts))
	#if info is not None:
		#gaVariant.info.
	#gaVariant.calls = rec.calls
	return gaVariant

def csMes(hdr):
	gaCallSet = variants_pb2.CallSet()
	gaCallSet.id = hdr.id
	gaCallSet.name = sample
	gaCallSet.bio_sample_id = sample


	return gaCallSet

#VCF file to read
file = sys.argv[1]
#file to output to
ofile = sys.argv[2]
sys.stdout = open(ofile, "w")
samp = pysam.VariantFile(file)
hdr = samp.header
#chromosome choice
chromo = "ref_brca1"
for rec in samp.fetch(chromo, 0, 100):
	print (json_format._MessageToJsonObject(vMes(rec), True))
for sample in hdr.samples:
	print (json_format._MessageToJsonObject(csmes(hdr), True))
