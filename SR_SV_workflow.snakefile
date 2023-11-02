from os.path import join
import pandas as pd
import glob

samples_bams = pd.read_table(config["sample_table"])
sample_bam_dict = { row.sample_name : row.bam_file for i,row in samples_bams.iterrows() }
#sample_cnv_dict = { row.sample_name : str(row.cnv_file)  for i,row in samples_bams.iterrows() }
sample_LRSV_dict = { row.sample_name : row.LR_SV for i,row in samples_bams.iterrows() if pd.notnull(row.LR_SV)}
LR_samples = list(sample_LRSV_dict.keys())
outdir = config["results_dir"]
samples = list(sample_bam_dict.keys())

localrules: all 

rule all:
	input:
		#expand("/oak/stanford/groups/euan/projects/tannerj/short_read_SV_workflow/results/merged/{sample}.hg38.diploidSV.merged.vcf.gz",
		#	sample = samples),
		#expand("/oak/stanford/groups/euan/projects/tannerj/short_read_SV_workflow/results/merged/{sample}.hg38.diploidSV.merged.stats.txt",
		#	sample = samples),
		#join(outdir, "merged", config["cohort_name"] + ".hg38.LR_overlap.counts.merged.txt"),
		join(outdir, config["cohort_name"] + ".COHORT_MERGED.hg38.short_read_wgs.SVs.vcf"),
		join(outdir, config["cohort_name"] + ".hg38.pop_SVs.paragraph_genotyped.vcf.gz"),

rule mantaSV:
	threads: 4
	resources:
		mem=500,
		time=24
	input:
		fasta = lambda w: config["references"][w.ref_name]['fasta'],
		fai = lambda w: config["references"][w.ref_name]['fasta']+ '.fai',
		bam = lambda w: sample_bam_dict[w.sample]
	params:
		tmp_dir = lambda w: join(config["project_dir"], "tmp/manta_" + w.sample + w.ref_name),
		scripts_dir = join(config["project_dir"], "helper_scripts"),
		out_dir = join(outdir, "mantaSV")
	output:
		vcf = join(outdir, "mantaSV/{sample}.{ref_name}.manta.diploidSV.vcf.gz"),
		tbi = join(outdir, "mantaSV/{sample}.{ref_name}.manta.diploidSV.vcf.gz.tbi")
	conda: 'envs/manta.yml' 
	shell: """
		ml samtools
		mkdir -p {params.tmp_dir}
		mkdir -p {params.out_dir}
		cd {params.tmp_dir}

		samtools view -H {input.bam} | grep '^@SQ' > bam_header.reference.txt 
		ln -s {input.fasta} genome.fa

		python {params.scripts_dir}/match_fasta_reference_to_bam.py \
			bam_header.reference.txt  \
			{input.fai} > genome.fa.fai

		# limit to larger chroms ( > 10MB)
		awk -v "OFS=\t" "\$2 > 10000000 || \$1 ~ /(M|MT)\$/ {{ print \$1,0,\$2 }}" genome.fa.fai | bgzip -c > cr.bed.gz

		tabix cr.bed.gz
		configManta.py --bam {input.bam} --referenceFasta genome.fa --runDir . --callRegions cr.bed.gz
		python2 ./runWorkflow.py -j {threads}
		python2 {params.scripts_dir}/convertInversion.py $(which samtools) {input.fasta} results/variants/diploidSV.vcf.gz \
			| bcftools view -O z -o {output.vcf}
		bcftools index --tbi {output.vcf}

		#rm -rf {params.tmp_dir}
	"""

rule dysgu_run: 
	threads: 4
	resources:
		mem=500,
		time=12
	input:
		fasta = lambda w: config["references"][w.ref_name]['fasta'],
		fai = lambda w: config["references"][w.ref_name]['fasta'] + '.fai',
		bam = lambda w: sample_bam_dict[w.sample]
	params:
		tmp_dir = lambda w: join("tmp/dysgu_" + w.sample + w.ref_name),
		max_depth = 200,
		out_dir = join(outdir, "dysgu")
	output:
		vcf = join(outdir, "dysgu/{sample}.{ref_name}.dysgu.diploidSV.vcf.gz"),
		tbi = join(outdir, "dysgu/{sample}.{ref_name}.dysgu.diploidSV.vcf.gz.tbi")
	conda: 'envs/dysgu.yml'	
	shell: """
		## clear and recreate working directory
		rm -rf {params.tmp_dir}
		mkdir -p {params.tmp_dir}
		mkdir -p {params.out_dir}
		cd {params.tmp_dir}

		## limit to larger chroms ( > 10MB)
		awk -v "OFS=\t" "\$2 > 10000000 || \$1 ~ /(M|MT)\$/ {{ print \$1,0,\$2 }}" {input.fai} > cr.bed

		dysgu run --clean \
				--pl pe --mode pe \
				-p {threads} \
				--max-cov {params.max_depth} \
				--thresholds 0.25,0.25,0.25,0.25,0.25 \
				--search cr.bed \
				{input.fasta} './temp' {input.bam} > {wildcards.sample}.tmp.vcf 
		bcftools sort -m 6G -O z -o {output.vcf} {wildcards.sample}.tmp.vcf
		bcftools index --tbi {output.vcf}

		#rm -rf {params.tmp_dir}
	"""
		

rule merge_within_sample:
	threads: 1
	resources:
		mem=6,
		time=2
	input:
		manta = rules.mantaSV.output[0],
		dysgu = rules.dysgu_run.output[0]
	params:
		out_dir = lambda w: join(outdir, 'merged'),
		tmp_dir = lambda w: join("tmp/merge_" + w.sample + w.ref_name),
		###cnv = lambda w: sample_cnv_dict[w.sample]
		cnv = 0,
		scripts_dir = join(config["project_dir"], "helper_scripts")
	output:
		vcf = join(outdir, "merged/{sample}.{ref_name}.diploidSV.merged.vcf"),
		stats = join(outdir, "merged/{sample}.{ref_name}.diploidSV.merged.stats.txt")
	conda: 'envs/manta.yml'
	shell: """
		
		mkdir -p {params.tmp_dir}

		cd {params.tmp_dir}
		export cnv_arg=""
		#if [ "{params.cnv}" != "nan" ]; then
		#	bcftools view -i 'QUAL>=20' -f PASS,cnvLength -U {params.cnv} > cnv.vcf
		#	export cnv_arg="cnv.vcf"
		#fi

		bcftools view -i 'QUAL>=60' -f PASS {input.manta} > manta.vcf
		bcftools view -f PASS {input.dysgu} > dysgu.vcf

		ls manta.vcf dysgu.vcf $cnv_arg > input.vcfs_list.txt
		SURVIVOR merge input.vcfs_list.txt 1000 1 1 1 0 30 tmp.vcf
		
		#variants (BND or TRA) that map to HLA contigs get confused in bcftools so we remove them
		grep -v 'HLA-' tmp.vcf | \
			python3 {params.scripts_dir}/consensus_genotype.py {wildcards.sample} | \
			bcftools sort -m 6G -Ov -o {output.vcf}
		bcftools index --tbi {output.vcf}
		
		SURVIVOR stats tmp.vcf -1 -1 -1 {output.stats}

		rm -rf {params.tmp_dir}
	"""


rule overlap_LR:
	threads: 1
	resources: 
		mem=12,
		time=2
	input:
		join(outdir, "merged/{sample}.{ref_name}.diploidSV.merged.vcf")
	params:
		LR = lambda w: sample_LRSV_dict[w.sample],
		region = config["stratify_region"],
	output:
		temp(join(outdir, "merged/{sample}.{ref_name}.LR_overlap.txt"))
	shell: """
		module load R/3.6.1
		Rscript helper_scripts/overlap.R {input} {params.LR} {params.region} {wildcards.sample} {output}
	"""
rule merge_overlaps:
	threads: 1
	resources:
		mem=4,
		time=1
	input: 
		expand(join(outdir, "merged/{sample}.{{ref_name}}.LR_overlap.txt"),
			sample=LR_samples)
	output:
		join(outdir, "merged", config["cohort_name"] + ".{ref_name}.LR_overlap.counts.merged.txt")
	shell: '''
		cat {input} | sort -k1,1r | uniq > {output}
	'''

rule SURVIVOR_pop_merge:
	threads: 1
	resources: 
		mem=24,
		time=6
	input: 
		expand(join(outdir, "merged/{sample}.{{ref_name}}.diploidSV.merged.vcf"), 
			sample = samples)
	output:
		vcf  = join(outdir, config["cohort_name"] + ".COHORT_MERGED.{ref_name}.short_read_wgs.SVs.vcf"),
		stats  = join(outdir, config["cohort_name"] + ".COHORT_MERGED.{ref_name}.short_read_wgs.SVs.stats")
	shell: '''
		ls {input} > tmp.vcf_list
		SURVIVOR merge tmp.vcf_list 1000 1 1 1 0 30 {output.vcf}
		SURVIVOR stats {output.vcf} -1 -1 -1 {output.stats}
		rm tmp.vcf_list
	'''
rule jasmine_population_merge:
	threads: 8
	resources:
		mem=18,
		time=8
	input:
		expand(join(outdir, "merged/{sample}.{{ref_name}}.diploidSV.merged.vcf"),
			sample = samples)
	params:
		tmp_dir = lambda w: join("tmp/jasmine_" + w.ref_name),
		scripts_dir = join(config["project_dir"], "helper_scripts")
	output:
		vcf = join(outdir, "jasmine", config["cohort_name"] + ".jasmine_merged.{ref_name}.vcf.gz"),
		tbi = join(outdir, "jasmine", config["cohort_name"] + ".jasmine_merged.{ref_name}.vcf.gz.tbi")
	conda: 'envs/jasmine_paragraph.yml'
	shell: """
		mkdir -p {params.tmp_dir}
		cd {params.tmp_dir} 

		echo "{input}" | tr ' ' '\n' > vcf_input.list
		jasmine -Xmx6g --threads {threads} --allow_intrasample --file_list=vcf_input.list  --out_file="jasmine.tmp.vcf"


		# add insertion seq, remove breakends (BND) which cant be genotyped, and sort
		bcftools filter -e'INFO/SVTYPE=="BND"||POS<=150||FILTER!="PASS"' jasmine.tmp.vcf | \
			python3 {params.scripts_dir}/add_insertion_seq.py | \
			bcftools sort -m 6G -Oz -o {output.vcf}

		#index
		bcftools index --tbi {output.vcf}

		rm -rf {params.tmp_dir}
	"""

rule paragraph_genotype:
	threads: 16
	resources:
		mem=18,
		time=8
	input:
		bam = lambda w: sample_bam_dict[w.sample],
		popvcf = rules.jasmine_population_merge.output[0],
		fasta = lambda w: config["references"][w.ref_name]
	params:
		tmp_dir = lambda w: join("tmp/paragraph_" + w.sample + w.ref_name),
		rd = config["mean_coverage"], # read depth
		rl = config["read_length"] #read length
	output:
		vcf = join(outdir, "paragraph/{sample}.{ref_name}.paragraph_genotyped.vcf.gz"),
		tbi = join(outdir, "paragraph/{sample}.{ref_name}.paragraph_genotyped.vcf.gz.tbi")
	conda: 'envs/jasmine_paragraph.yml'
	shell: """
		mkdir -p {params.tmp_dir}
		cd {params.tmp_dir} 


		# limit to main chroms (>10Mb)
		awk -v "OFS=\t" "\$2 > 10000000 || \$1 ~ /(M|MT)\$/ {{ print \$1,0,\$2 }}" {input.fasta}.fai > cr.bed
	   	bedtools intersect -header -a {input.popvcf} -b cr.bed > input.cr.vcf

		echo -e "id\tpath\tdepth\tread length" > sample.manifest
		echo -e "{wildcards.sample}\t{input.bam}\t{params.rd}\t{params.rl}" >> sample.manifest
		M=$(({params.rd} * 6))

		# this is the main paragraph entrypoint
		multigrmpy.py -i input.cr.vcf \
			-m sample.manifest \
			-r {input.fasta} \
			-o tmp/ \
			-t {threads} \
			-M $M

		# duphold adds depth annotations looking at coverage fold-change around Svs
		#duphold -d -v tmp/genotypes.vcf.gz -b {input.bam} -f {input.fasta} -t {threads} -o {output.vcf}
		bcftools sort -Oz -o {output.vcf} tmp/genotypes.vcf.gz
		bcftools index --tbi {output.vcf}

		rm -rf {params.tmp_dir}
	"""
		
rule merge_genotypes:
	threads: 1
	resources:
		mem=12,
		time=4
	input:
		expand(join(outdir, "paragraph/{sample}.{{ref_name}}.paragraph_genotyped.vcf.gz"),
			sample = samples)
	output:
		vcf = join(outdir, config["cohort_name"] + ".{ref_name}.pop_SVs.paragraph_genotyped.vcf.gz"),
		tbi = join(outdir, config["cohort_name"] + ".{ref_name}.pop_SVs.paragraph_genotyped.vcf.gz.tbi")
	conda: 'envs/jasmine_paragraph.yml'
	shell: """
		bcftools merge -Oz -o {output.vcf} {input}
		bcftools index -f --tbi {output.vcf}
	"""
