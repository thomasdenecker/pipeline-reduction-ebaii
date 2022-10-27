SAMPLES, = glob_wildcards(config["dataDir"]+"{sample}_R1.fastq.gz")
IDX = ["1","2","3","4","5","6","7","8"]

rule all :
  input :
    expand("FastQC/{sample}_R1_fastqc.html", sample=SAMPLES),
    expand("FastQC/{sample}_R2_fastqc.html", sample=SAMPLES),
    expand(config["dataRes"]+"{sample}_R{R}_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq.gz", sample=SAMPLES, R=["1","2"])


rule compress:
  output:
    config["dataRes"]+"{sampleRed}"+config["chr"]+"_"+str(config["nbAlign"])+".fastq.gz"
  input:
    config["dataRes"]+"{sampleRed}"+config["chr"]+"_"+str(config["nbAlign"])+".fastq"
  log:
    "Logs/compress_{sampleRed}"+config["chr"]+"_"+str(config["nbAlign"])+".log"
  shell:
    "gzip {input}"    


rule create_fastq:
  output:
    R1=config["dataRes"]+"{sample}_R1_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq",
    R2=config["dataRes"]+"{sample}_R2_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq" 
  input:
    "Tmp/{sample}_select.sam"
  log:
    "Logs/{sample}_createFastqRed.log"
  conda:
    "samtool.yaml"
  shell:
    "samtools fastq -n -1 {output.R1} -2 {output.R2} Tmp/{wildcards.sample}_select.sam 2> {log} "


rule extract_reads:
  output:
    "Tmp/{sample}_select.sam",
  input:
    "Tmp/{sample}.bam"
  log:
    "Logs/{sample}_extract.log"
  conda:
    "samtool.yaml"
  shell:
    "samtools view -H {input} > {output} 2> {log} ; "
    "set +o pipefail ; "
    "samtools view {input} | head -n "+str(config["nbAlign"])+" | grep -v "+config["chr"]+" >> {output} 2>> {log} ; "
    "samtools view {input} | grep "+config["chr"]+" >> {output} 2>> {log} ; "


rule hisat2_mapping:
  output:
    temp("Tmp/{sample}.bam")
  input:
#    R1=config["dataDir"]+"{sample}.fastq.gz" if config["rnaType"]=="single-end" else ??
    expand("Tmp/GenomeIdx/GenomeIdx.{idx}.ht2", idx=IDX),
    R1="Tmp/{sample}_R1.fastq",
    R2="Tmp/{sample}_R2.fastq"
  threads: config["threads_mapping"]
  resources:
    cpus=config["threads_mapping"]
  log:
    "Logs/{sample}_bwt2_mapping.log"
  conda: 
    "mapping.yml"
  shell:
    "hisat2 -x Tmp/GenomeIdx/GenomeIdx -p {threads} -k 1 --no-mixed --rna-strandness FR -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -@ {threads} -b -o {output}"
#  run: 
#    """
#     if (config["rnaType"]=="paired-end"):
#        shell("hisat2 -x Tmp/GenomeIdx -1 {input.R1} -2 {input.R2} -S {output} 2> {log}")
#     else:
#        shell("hisat2 -x Tmp/GenomeIdx -U {input.R1} -S {output} 2> {log}")
#    """
 

rule genome_hisat2_index:
  output:
    expand("Tmp/GenomeIdx/GenomeIdx.{ext}.ht2", ext=IDX)
  input:
    fna=config["genome"]
  params:
    idx="Tmp/GenomeIdx/GenomeIdx"
  resources:
    mem_mb=config["mem_index"]
  conda:
    "hisat2.yml"
  log:
    log1="Logs/genome_hisat2_index.log1",
    log2="Logs/genome_hisat2_index.log2"
  shell: 
    "hisat2-build {input.fna} {params.idx} 1>{log.log1} 2>{log.log2}"


rule fastqc:
  output: 
    "FastQC/{sample}_fastqc.zip",
    "FastQC/{sample}_fastqc.html"
  input: 
    "Tmp/{sample}.fastq"
  log:
    log1="Logs/{sample}_fastqc.log1",
    log2="Logs/{sample}_fastqc.log2"
  conda: "condaEnv4SmkRules/fastqc.yml"
  shell: "fastqc --outdir FastQC/ {input} 1>{log.log1} 2>{log.log2}"


rule uncompress:
  output:
    temp("Tmp/{sample}.fastq")
  input:
    config["dataDir"]+"{sample}.fastq.gz"
  log:
    "Logs/{sample}_uncompress.log"
  run:
    shell("gunzip -c {input} > {output} ")



