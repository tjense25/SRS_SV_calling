# SRS_SV_calling
snakemake workflow to call structural variants from short-read whole genome sequencing data

# Short Read Structural Variant Calling
This workflow uses Manta and Dysgu to call structural variants from shrot read bams, then merges these calls with jasmine to create a population-level reference of SVs, and finally uses paragraph to genotype in each sample the population reference SVs, allowing for accurate genotypes and robust allele frequency estimation. 

# Dependencies 
This workflow requires, conda, snakemake, and mamba as dependencies. Snakemake and mamba can be installed using conda with the following commands. 
```
conda install -c conda-forge mamba 
conda install -c bioconda snakemake
```

# To Run Workflow:
1. clone this git repository into a working directory
```
git clone https://github.com/tjense25/SRS_SV_calling.git
```

2. Update SR_SV_workflow.config.yml to point to the correct data

* **cohort_name**: name of cohort used to label aggregated SV vcf file

* **read_length**: readlength of short-reads (usually 150)

* **read_depth**: approximate average depth of coverage of short-read bams. need not be precise, err on the side of being too high

* **sample_table**: path to a tab-delimited, two column file containing sample name and short-read bam path for all samples (example sample_table file - **sample_bams.txt**)

* **project_dir**: full path of the working directory that contains the snakefile

* **results_dir**: full path to the directory where results should be written

* **references**: **hg38**: full path to hg38 reference genome .fasta file. Must be faidx indexed. Must also match reference used for short-read alignment exactly. 

* **threads**: number of thread to use for multi-threading parallelization 

3. To run snakemake on Sherlock cluster using SLURM job scheduler usee the following command:
```
snakemake --profile sherlock --jobs 100 --snakefile SR_SV_workflow.snakefile --configfile SR_SV_workflow.config.yml --use-conda 
```

To run snakemake locally on a server:
```
snakemake --jobs 16 --snakefile SR_SV_workflow.snakefile --configfile SR_SV_workflow.config.yml --use-conda
```
