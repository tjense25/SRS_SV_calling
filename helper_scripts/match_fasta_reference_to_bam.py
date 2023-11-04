import sys

bam_header_file = sys.argv[1]
reference_index = sys.argv[2]
seqs_in_bam = set()

with open(bam_header_file, 'r') as in_file:
	for line in in_file:
		toks = line.strip().split('\t')
		seqs_in_bam.add(toks[1].split('SN:')[1])

with open(reference_index, 'r') as in_file:
	for line in in_file:
		toks = line.strip().split('\t')
		seq = toks[0]
		if seq in seqs_in_bam:
			print(line.strip())
