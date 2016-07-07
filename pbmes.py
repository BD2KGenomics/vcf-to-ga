import variants_pb2
import sys
import pysam
from pysam import VariantFile
import google.protobuf.json_format as json_format




def varmes(rec):
	gaVariant = variants_pb2.Variant()
	gaVariant.id = rec.id
	gaVariant.reference_name = rec.contig
	#gaVariant.variant_set_id = rec.variant_set_id
	#gaVariant.names
	#gaVariant.created = rec.created
	#gaVariant.updated = rec.updated
	gaVariant.start = rec.start
	gaVariant.end = rec.stop
	gaVariant.reference_bases = rec.ref
	if rec.alts is not None:
		gaVariant.alternate_bases.extend(list(rec.alts))
	#if rec.info is not None:
		#gaVariant.info.append(list(rec.info))
	#gaVariant.calls = rec.calls
	return gaVariant

file = sys.argv[1]
ofile = sys.argv[2]
sys.stdout = open(ofile, "w")
samp = pysam.VariantFile(file)
chrs = "ref_brca1"
#start = int(raw_input("start: "))
#stop = int(raw_input("stop: "))
for rec in samp.fetch(chrs, 0, 100):
	print (json_format._MessageToJsonObject(varmes(rec), True))

