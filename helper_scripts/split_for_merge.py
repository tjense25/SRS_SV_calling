#!/bin/bash/python
import sys
import os

infile = sys.argv[1]
outdir=sys.argv[2]
outfile=os.path.join(outdir, os.path.basename(infile))

chrom_map = { 
		"chr1" : outfile.replace(".vcf", ".chr1.vcf"),
		"chr2" : outfile.replace(".vcf", ".chr2.vcf"),
		"chr3" : outfile.replace(".vcf", ".chr3.vcf"),
		"chr4" : outfile.replace(".vcf", ".chr4.vcf"),
		"chr5" : outfile.replace(".vcf", ".chr5.vcf"),
		"chr6" : outfile.replace(".vcf", ".chr6.vcf"),
		"chr7" : outfile.replace(".vcf", ".chr7.vcf"),
		"chr8" : outfile.replace(".vcf", ".chr8.vcf"),
		"chr9" : outfile.replace(".vcf", ".chr9.vcf"),
		"chr10" :outfile.replace(".vcf",  ".chr10.vcf"),
		"chr11" :outfile.replace(".vcf",  ".chr11.vcf"),
		"chr12" :outfile.replace(".vcf",  ".chr12.vcf"),
		"chr13" :outfile.replace(".vcf",  ".chr13.vcf"),
		"chr14" :outfile.replace(".vcf",  ".chr14.vcf"),
		"chr15" :outfile.replace(".vcf",  ".chr15.vcf"),
		"chr16" :outfile.replace(".vcf",  ".chr16.vcf"),
		"chr17" :outfile.replace(".vcf",  ".chr17.vcf"),
		"chr18" :outfile.replace(".vcf",  ".chr18.vcf"),
		"chr19" :outfile.replace(".vcf",  ".chr19.vcf"),
		"chr20" :outfile.replace(".vcf",  ".chr20.vcf"),
		"chr21" :outfile.replace(".vcf",  ".chr21.vcf"),
		"chr22" :outfile.replace(".vcf",  ".chr22.vcf"),
		"chrX" : outfile.replace(".vcf", ".chrX.vcf"),
		"chrY" : outfile.replace(".vcf", ".chrY.vcf")}
EXTRA = outfile.replace(".vcf", ".chrEXTRA.vcf")

for f in chrom_map:
	chrom_map[f] = open(chrom_map[f], 'w')
EXTRA = open(EXTRA, 'w')

infile = open(infile,'r')

for line in infile:
	if line.startswith('#'):
		for f in chrom_map:
			print(line.strip(), file=chrom_map[f])
		print(line.strip(), file=EXTRA)
	COLS=line.strip().split("\t")
	if COLS[0] in chrom_map:
		print(line.strip(),file = chrom_map[COLS[0]])
	else:
		print(line.strip(),file=EXTRA)

infile.close()
for f in chrom_map:
	chrom_map[f].close()
EXTRA.close()
