# vcf_to_ga
converts VCF into GA4GH protobuf messages

To begin:
Turn on a virtual environment and pip install requirements

python scr.py -i vcfFile -db DatabaseName

When running the program, input file (-i filename) is required, specifying database (-db databasenameyouwant) is optional.  If you do not specify -db, only directory structured output will be created.  You can also specify chromosome with -ch chromosomeName, this is optional and if you do not specify, the entire file will be fetched.
