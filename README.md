# vcf_to_ga
converts VCF into GA4GH

To begin:
Turn on a virtual environment and pip install requirements

python scr.py -i vcfFile 

optional commands:
-db DatabaseName
-ch chromosomeName
-pb AnyTextHereWorks

When running the program, input file (-i filename) is required, specifying database (-db databasenameyouwant) is optional.  If you do not specify -db, only directory structured output will be created.  You can also specify chromosome with -ch chromosomeName, this is optional and if you do not specify, the entire file will be fetched.  If you wish to have .pb as output instead of .txt into the directory use the extension -pb pb.  If this is specified it will output only to .pb.  If this is not specified it will output to .txt.  This command is optional and can be used with or without outputting to database.
