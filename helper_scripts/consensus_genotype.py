#!/usr/bin/python
from collections import defaultdict
import sys


sample_name = sys.argv[1] 
for line in sys.stdin:
	if line.startswith('##'): 
		print(line.strip())
		continue

	toks = line.strip().split()

	if line.startswith('#'): #Add consensus column as first sample column
		toks.insert(9, sample_name + '_shortread')
		print('\t'.join(toks))
		continue 

	gts = [ x.split(':')[0] for x in toks[9:len(toks)] ]
	
	if gts.count('0/1') >= 2: #if majority of genotypes are 0/1 put the first 0/1 GT field as consensus field
		for i in range(len(gts)):
			if gts[i] == '0/1':
				toks.insert(9, toks[9 + i])
				break
	elif gts.count('1/1') >= 2: # if majority of genotypes are 1/1 put the first 1/1 GT field as consensus field
		for i in range(len(gts)):
			if gts[i] == '1/1':
				toks.insert(9, toks[9 + i])
				break
	elif gts.count('0/1') == 1 and gts.count('1/1') == 1: #if one GT is 0/1 and another is 1/1, make the 0/1 genotype the consensus GT
		for i in range(len(gts)):
			if gts[i] == '0/1':
				toks.insert(9, toks[9 + i])
				break

	else: #in there is only one 0/1 or 1/1 genotype (rest our ./. or 0/0) so we do not include this variant in consensus vcf
		toks.insert(9, toks[9].replace('./.', '0/0'))

	print('\t'.join(toks))

