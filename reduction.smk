SAMPLES, = glob_wildcards(config["dataDir"]+"{sample}_R1.fastq.gz")
IDX = ["1","2","3","4","5","6","7","8"]

rule all :
  input :
    expand("FastQC/{sample}_R1_fastqc.html", sample=SAMPLES),
    expand("FastQC/{sample}_R2_fastqc.html", sample=SAMPLES),
    expand(config["dataDir"]+"{sample}_R{R}_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq.gz", sample=SAMPLES, R=["1","2"])


rule compress:
  output:
    "{sample}.fastq.gz"
  input:
    "{sample}.fastq"
  log:
    "Logs/{sample}_compress.log"
  shell:
    "gzip {input}"    


rule uncompress:
  output:
    "{sample}.fastq"
  input:
    "{sample}.fastq.gz"
  log:
    "Logs/{sample}_uncompress.log"
  shell:
    "gunzip {input}"    


rule extract_reads:
  output:
    config["dataDir"]+"{sample}_R1_red.fastq",
    config["dataDir"]+"{sample}_R2_red.fastq"
  input:
    "Tmp/{sample}.bam"
  log:
    "Logs/{sample}_extract.log"
  conda:
    "samtool.yaml"
  shell:
    "samtools view -H {input} > Tmp/{sample}_select.sam "
    "samtools view {input} | head -n "+str(config["nbAlign"])+" | grep -v "+config["chr"]+" >> Tmp/{sample}_select.sam "
    "samtools view {input} | grep "+config["chr"]+" >> Tmp/{sample}_select.sam "
    "samtools fastq -n -1 "+config["dataDir"]+"{sample}_R1_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq -2 "+config["dataDir"]+"{sample}_R2_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq {input} "


rule hisat2_mapping:
  output:
    "Tmp/{sample}.bam"
  input:
#    R1=config["dataDir"]+"{sample}.fastq.gz" if config["rnaType"]=="single-end" else ??
    expand("Tmp/GenomeIdx.{idx}.ht2", idx=IDX),
    R1=config["dataDir"]+"{sample}_R1.fastq",
    R2=config["dataDir"]+"{sample}_R2.fastq"
  log:
    "Logs/{sample}_bwt2_mapping.log"
  conda: 
    "mapping.yml"
  shell:
    "hisat2 -x Tmp/GenomeIdx -k 1 --no-mixed --rna-strandness FR -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -b -o {output}"
#  run: 
#    """
#     if (config["rnaType"]=="paired-end"):
#        shell("hisat2 -x Tmp/GenomeIdx -1 {input.R1} -2 {input.R2} -S {output} 2> {log}")
#     else:
#        shell("hisat2 -x Tmp/GenomeIdx -U {input.R1} -S {output} 2> {log}")
#    """
 

rule genome_hisat2_index:
  output:
    expand("Tmp/GenomeIdx.{ext}.bt2", ext=IDX)
  input:
    fna=config["dataDir"]+config["genome"]+".fna",
    idx=expand("Tmp/GenomeIdx.{idx}.ht2", idx=IDX)
  conda:
    "hisat2.yml"
  log:
    log1="Logs/genome_hisat2_index.log1",
    log2="Logs/genome_hisat2_index.log2"
  shell: 
    "hisat2-build {input.fna} {input.idx} 1>{log.log1} 2>{log.log2}"

rule fastqc:
  output: 
    "FastQC/{sample}_fastqc.zip",
    "FastQC/{sample}_fastqc.html"
  input: 
    config["dataDir"]+"{sample}.fastq.gz"
  log:
    log1="Logs/{sample}_fastqc.log1",
    log2="Logs/{sample}_fastqc.log2"
  conda: "condaEnv4SmkRules/fastqc.yml"
  shell: "fastqc --outdir FastQC/ {input} 1>{log.log1} 2>{log.log2}"

