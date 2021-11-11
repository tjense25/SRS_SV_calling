import sys
import signal

# Don't turn these signal into exceptions, just die.
signal.signal(signal.SIGINT, signal.SIG_DFL)
signal.signal(signal.SIGPIPE, signal.SIG_DFL)

header={}
for line in sys.stdin:
	if line.startswith("##"): 
		sys.stdout.write(line)
	elif line.startswith("#CHROM"):
		header = { name: i for i,name in enumerate(line[1:].strip().split("\t")) }
		sys.stdout.write(line)
	else:
		cols = line.strip().split("\t")
		info = { key_val[0] : None if len(key_val) == 1 else key_val[1] for key_val in map(lambda x: x.split('='), cols[header["INFO"]].split(';'))}
		if cols[header["REF"]] == '.': cols[header["REF"]] = "N"
		if cols[header["ALT"]] == "<INS>":
			sv_len = int(info["SVLEN"])
			left = info['LEFT_SVINSSEQ'] if 'LEFT_SVINSSEQ' in info else ''
			right = info['RIGHT_SVINSSEQ'] if 'RIGHT_SVINSSEQ' in info else ''
			middle = info['EXPSEQ'] if 'EXPSEQ' in info and info['EXPSEQ'] else 'N'
			ins_seq = cols[header["REF"]][0] + left + middle*((sv_len - len(left) - len(right)) // len(middle)) + right
			cols[header["ALT"]] = ins_seq.upper()
		elif cols[header["ALT"]] == "<INV>":
			if int(info["SVLEN"]) == 0:
				cols[header["ALT"]] = info["SVINSSEQ"]
				info["SVTYPE"] = "INS"
				info["SVLEN"]=str(len(info["SVINSSEQ"]))
				cols[header["INFO"]] = ";".join([ k if not v else k+'='+v for k,v in info.items()])



		print('\t'.join(cols), file=sys.stdout, flush=True)
