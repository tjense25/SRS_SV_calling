from os.path import join
import pandas as pd
import glob

samples_bams = pd.read_table(config["sample_table"])
sample_bam_dict = { row.sample_name : row.bam_file for i,row in samples_bams.iterrows() }
outdir = config["results_dir"]

localrules: all 

samples=['UDN467147','UDN054366']
rule all:
	input:
		join(outdir, config["cohort_name"] + ".hg38.pop_SVs.paragraph_genotyped.vcf.gz"),

rule mantaSV:
	threads: 1
	resources:
		mem=32,
		time=24
	input:
		fasta = lambda w: config["references"][w.ref_name],
		fai = lambda w: config["references"][w.ref_name] + '.fai',
		bam = lambda w: sample_bam_dict[w.sample]
	params:
		tmp_dir = lambda w: join(config["project_dir"], "tmp/manta_" + w.sample + w.ref_name),
		scripts_dir = join(config["project_dir"], "helper_scripts")
	output:
		vcf = join(outdir, "mantaSV/{sample}.{ref_name}.manta.diploidSV.vcf.gz"),
		tbi = join(outdir, "mantaSV/{sample}.{ref_name}.manta.diploidSV.vcf.gz.tbi")
	conda: 'envs/manta.yml' 
	shell: """
		mkdir -p {params.tmp_dir}
		cd {params.tmp_dir}

		# limit to larger chroms ( > 10MB)
		awk -v "OFS=\t" "\$2 > 10000000 || \$1 ~ /(M|MT)\$/ {{ print \$1,0,\$2 }}" {input.fai} | bgzip -c > cr.bed.gz

		tabix cr.bed.gz
		configManta.py --bam {input.bam} --referenceFasta {input.fasta} --runDir . --callRegions cr.bed.gz
		python2 ./runWorkflow.py -j {threads}
		python2 {params.scripts_dir}/convertInversion.py $(which samtools) {input.fasta} results/variants/diploidSV.vcf.gz \
			| bcftools view -O z -o {output.vcf}
		bcftools index --tbi {output.vcf}

		rm -rf {params.tmp_dir}
	"""

rule dysgu_run: 
	threads: 1
	resources:
		mem=32,
		time=24 
	input:
		fasta = lambda w: config["references"][w.ref_name],
		fai = lambda w: config["references"][w.ref_name] + '.fai',
		bam = lambda w: sample_bam_dict[w.sample]
	params:
		tmp_dir = lambda w: join("tmp/dysgu_" + w.sample + w.ref_name),
		max_depth = 200
	output:
		vcf = join(outdir, "dysgu/{sample}.{ref_name}.dysgu.diploidSV.vcf.gz"),
		tbi = join(outdir, "dysgu/{sample}.{ref_name}.dysgu.diploidSV.vcf.gz.tbi")
	conda: 'envs/dysgu.yml'	
	shell: """
		## clear and recreate working directory
		rm -rf {params.tmp_dir}
		mkdir -p {params.tmp_dir}
		cd {params.tmp_dir}

		# limit to larger chroms ( > 10MB)
		awk -v "OFS=\t" "\$2 > 10000000 || \$1 ~ /(M|MT)\$/ {{ print \$1,0,\$2 }}" {input.fai} > cr.bed

		dysgu run --clean \
				--pl pe --mode pe \
				-p {threads} \
				--max-cov {params.max_depth} \
				--thresholds 0.25,0.25,0.25,0.25,0.25 \
				--search cr.bed \
				{input.fasta} './temp' {input.bam} > {wildcards.sample}.tmp.vcf 
		bcftools view -O u -o {wildcards.sample}.R.bcf  {wildcards.sample}.tmp.vcf
		bcftools sort -m 6G -O z -o {output.vcf} {wildcards.sample}.R.bcf
		bcftools index --tbi {output.vcf}

		rm -rf {params.tmp_dir}
	"""
		

rule merge_within_sample:
	threads: 1
	resources:
		mem=6,
		time=2
	input:
		rules.mantaSV.output[0],
		rules.dysgu_run.output[0]
	output:
		join(outdir, "merged/{sample}.{ref_name}.diploidSV.merged.vcf")
	conda: 'envs/manta.yml'
	shell: """
		 bcftools concat -a -O v -o {output} {input}
	"""

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

		#rm -rf {params.tmp_dir}
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


