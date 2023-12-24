from os.path import join
import pandas as pd
import glob

localrules: all
wildcard_constraints:
  source="consensus_SVs.jasmine_irisRefined|SVs.sniffles|sniffles_jointcall|spectre_CNV"

               
sample_table={}
sample_vcf={}
sample_snf={}
cohort_size = {}
for ref_name in config["references"]:
    sample_table[ref_name] = pd.read_table(config["references"][ref_name])
    sample_vcf[ref_name] = {"consensus_SVs.jasmine_irisRefined": {},
                            "SVs.sniffles": {}, 
                            "spectre_CNV":{}}
    cohort_size[ref_name] = 0
    sample_snf[ref_name] = {}
    for i,row in sample_table[ref_name].iterrows():
        sample_vcf[ref_name]["consensus_SVs.jasmine_irisRefined"][row.SampleName] = row.Consensus_VCF
        sample_vcf[ref_name]["SVs.sniffles"][row.SampleName] = row.sniffles_VCF
        sample_vcf[ref_name]["spectre_CNV"][row.SampleName] = row.Spectre_CNV
        sample_snf[ref_name][row.SampleName] = row.sniffles_SNF
        cohort_size[ref_name] += 1

scratch="/tmp/tannerj"
rule all:      
    input:     
        #expand("ALNG.{ref_name}.{source}.AF_annotated.vcf.gz",
        #       ref_name = "GRCh38",
        #       source = ["SVs.sniffles", "spectre_CNV", "consensus_SVs.jasmine_irisRefined"]),
        #       #source = ["consensus_SVs.jasmine_irisRefined", "SVs.sniffles"]),
        #"ALNG.GRCh38.sniffles_jointcall.AF_annotated.vcf.gz"
        "ALNG.GRCh38.consensus_SVs.jasmine_irisRefined.AF_annotated.vcf"

#chroms = ["chr" + str(x) for x in list(range(1,23)) + ["X","Y", "EXTRA"]]
chroms = ["chr" + str(x) for x in list(range(1,23)) + ["X","Y"]]

rule sniffles_joint_genotype:
    threads: 16
    resources:
        mem=256,
        time=24
    input:
        lambda w: sample_snf[w.ref_name].values()
    params:
        snf_tsv = "tmp/ALNG.snf.tsv",
        tmp_vcf = "tmp/ALNG.{ref_name}.sniffles_jointcall.vcf"
    output:
        vcf = "ALNG.{ref_name}.sniffles_jointcall.vcf.gz",
    conda: "envs/LR_workflow.yaml"
    shell: """
        ls {input} > {params.snf_tsv}
        sniffles -t 16 --input {params.snf_tsv} --vcf {params.tmp_vcf}
        bcftools sort -Oz -o {output.vcf} {params.tmp_vcf}
        rm {params}
    """

rule split_chroms_for_merge:
    threads: 1 
    resources:
        mem=12,
        time=4
    input:
        lambda w: sample_vcf[w.ref_name][w.source][w.sample] 
    params:
        outdir="tmp"
    output:
        temp(expand("tmp/{{sample}}.{{ref_name}}.{{source}}.{chrom}.vcf",
            chrom=chroms))
    shell: """
        python helper_scripts/split_for_merge.py {input} {params.outdir}
    """

rule jasmine_merge_all:
    threads: 16
    resources:
        mem=112,
        time=32
    input:
        vcfs = lambda w: expand("tmp/{sample}.{{ref_name}}.{{source}}.{{chrom}}.vcf",
            sample=sample_vcf[w.ref_name][w.source].keys()),
        ref = lambda w: config["reference_params"][w.ref_name]["fasta"]
    params:
        tmp_dir = join(scratch, "jasmine_{chrom}_{ref_name}_{source}"),
        vcf_list = join(scratch, "jasmine_{chrom}_{ref_name}_{source}/SV_vcfs.{ref_name}.{chrom}.list")
    output:
        temp("tmp/ALNG.{ref_name}.{source}.{chrom}.vcf")
    conda: "envs/jasmine.yaml"
    shell: """
        mkdir -p {params.tmp_dir}
        module load java/18.0.2.1
        ls {input.vcfs} > {params.vcf_list}
        jasmine -Xmx96g -Xms24g --ignore_strand --file_list={params.vcf_list} \
                        --normalize_type --dup_to_ins \
                        --default_zero_genotype --centroid_merging \
                        out_dir={params.tmp_dir} genome_file={input.ref} --output_genotypes \
                        out_file={output} threads={threads}
        rm -rf {params.tmp_dir}
    """

rule gather_jasmine:
    threads: 1
    resources: 
        mem=24,
        time=8
    input:
        expand("tmp/ALNG.{{ref_name}}.{{source}}.{chrom}.vcf",
            chrom=chroms)
    params:
        tmp_vcf="tmp/ALNG.cohort_combined.{ref_name}.{source}.vcf"
    output:
        vcf="ALNG.{ref_name}.{source}_jasmine.vcf.gz"
    conda: "envs/LR_workflow.yaml"
    shell: """
        bcftools concat --ligate-force -o {params.tmp_vcf} {input}
        cat {params.tmp_vcf} | python helper_scripts/prep_for_sort.py | bcftools sort -Oz -o {output.vcf} -
    """

rule annotate_AFs:
    threads: 1
    resources:
        mem=96,
        time=24
    input:
        "ALNG.{ref_name}.{source}.vcf.gz"
    output:
        "ALNG.{ref_name}.{source}.AF_annotated.txt"
    params: 
      cohort_size=lambda w: cohort_size[w.ref_name]
    conda: "r"
    shell: """
        Rscript annotate_AFs.R {input} {output} {params.cohort_size}
    """

rule convert_to_vcf_unzip:
    threads: 1
    resources:
        mem=24,
        time=12
    input:
        "{prefix}.AF_annotated.txt",
        "{prefix}.vcf.gz"
    output:
        vcf = "{prefix}.AF_annotated.vcf",
    conda: "envs/LR_workflow.yaml"
    shell: """
        python transfer_afs.py {input} > {output.vcf}
    """

rule convert_to_vcf:
    threads: 1
    resources:
        mem=24,
        time=12
    input:
        "{prefix}.AF_annotated.txt",
        "{prefix}.vcf.gz"
    output:
        vcf = "{prefix}.AF_annotated.vcf.gz",
        tbi = "{prefix}.AF_annotated.vcf.gz.tbi"
    conda: "envs/LR_workflow.yaml"
    shell: """
        python transfer_afs.py {input} | \
            bgzip -c > {output.vcf}
        bcftools index --tbi {output.vcf}
    """

rule merge_spectre:
  threads: 1
  resources:
    mem=24,
    time=12,
  input:
    lambda w: sample_vcf[w.ref_name]["spectre_CNV"].values()
  params:
    vcf_list = lambda w: "ALNG." + w.ref_name + ".spectre_CNV.list",
    tmp_output = "ALNG.{ref_name}.spectre_CNV.vcf" 
  output:
    "ALNG.{ref_name}.spectre_CNV.vcf.gz"
  conda: "envs/LR_workflow.yaml"
  shell: """
    ls {input} > {params.vcf_list}
    SURVIVOR merge {params.vcf_list} 1000 1 1 1 0 10000 {params.tmp_output}
    bcftools sort -Oz -o {output} {params.tmp_output}
    rm {params.tmp_output}
  """
