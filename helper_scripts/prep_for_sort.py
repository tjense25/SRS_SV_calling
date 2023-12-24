#!/usr/bin/python
from collections import defaultdict
import sys


for line in sys.stdin:
	if line.startswith('##'): 
		print(line.strip())
		continue

	COLS = line.strip().split("\t")
	if line.startswith("#CHROM"):
		HEADER = { name: i for i,name in enumerate(line[1:].strip().split("\t")) }
		print('\t'.join(COLS))
		continue

	INFO = { key_val[0] : None if len(key_val) == 1 else key_val[1] for key_val in map(lambda x: x.split('='), COLS[HEADER["INFO"]].split(';'))}

	INFO_OUTFIELDS=["SVTYPE","SVLEN","END", "CHR2","END2","STRAND","READ_SUPPORT","SVMETHOD","IRIS_PROCESSED","IRIS_REFINED", "SUPP", "SUPP_VEC"]
	INFO_OUT ="PRECISE;" if "PRECISE" in INFO else "IMPRECISE;"
	INFO_OUT += ";".join(['='.join([x,INFO[x]]) for x in INFO_OUTFIELDS if x in INFO])
	COLS[HEADER["INFO"]]=INFO_OUT
	print('\t'.join(COLS))
