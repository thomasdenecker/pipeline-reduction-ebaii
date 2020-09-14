SAMPLES, = glob_wildcards(config["dataDir"]+"{sample}_R1.fastq.gz")
IDX = ["1","2","3","4","5","6","7","8"]

rule all :
  input :
    expand("FastQC/{sample}_R1_fastqc.html", sample=SAMPLES),
    expand("FastQC/{sample}_R2_fastqc.html", sample=SAMPLES),
    expand(config["dataDir"]+"{sample}_R{R}_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq", sample=SAMPLES, R=["1","2"])


# rule compress:
#  output:
#    expand(config["dataDir"]+"{sample}_R{R}_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq.gz", sample=SAMPLES, R=["1","2"])
#  input:
#    expand(config["dataDir"]+"{sample}_R{R}_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq", sample=SAMPLES, R=["1","2"])
#  log:
#    "Logs/compress.log"
#  shell:
#    "gzip {input}"    


rule uncompress:
  output:
    "Tmp/{sample}.fastq"
  input:
    config["dataDir"]+"{sample}.fastq.gz"
  log:
    "Logs/{sample}_uncompress.log"
  run:
    shell("gunzip -c {input} > {output} ")
#rule uncompress:
#  output:
#    R1="Tmp/{sample}_R1.fastq",
#    R2="Tmp/{sample}_R2.fastq"
#  input:
#    R1=config["dataDir"]+"{sample}_R1.fastq.gz",
#    R2=config["dataDir"]+"{sample}_R2.fastq.gz"
#  log:
#    "Logs/{sample}_uncompress.log"
#  run:
#    shell("gunzip -c {input.R1} > {output.R1} ")
#    shell("gunzip -c {input.R2} > {output.R2} ")

rule create_fastq:
  output:
    R1=config["dataDir"]+"{sample}_R1_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq",
    R2=config["dataDir"]+"{sample}_R2_"+config["chr"]+"_"+str(config["nbAlign"])+".fastq" 
  input:
    "Tmp/{sample}_selectok.sam"
  log:
    "Logs/{sample}_extract.log"
  conda:
    "samtool.yaml"
  shell:
    "samtools fastq -n -1 {output.R1} -2 {output.R2} Tmp/{wildcards.sample}_select.sam "

rule extract_reads3:
  output:
    select="Tmp/{sample}_select3.sam",
    end="Tmp/{sample}_selectok.sam"
  input:
    "Tmp/{sample}_select2ok.sam"
  log:
    "Logs/{sample}_extract.log"
  conda:
    "samtool.yaml"
  shell:
    "samtools view {input} | grep "+config["chr"]+" >> {output.select} ; "
    "cp {output.select} {output.end} "

rule extract_reads2:
  output:
    select="Tmp/{sample}_select2.sam",
    end="Tmp/{sample}_select2ok.sam"
  input:
    "Tmp/{sample}.bam"
  log:
    "Logs/{sample}_extract.log"
  conda:
    "samtool.yaml"
  shell:
    "set +e; samtools view {input} | head -n "+str(config["nbAlign"])+" | grep -v "+config["chr"]+" >> {output.select} ; "
    "cp {output.select} {output.end} "


rule extract_reads1:
  output:
    select="Tmp/{sample}_select1.sam",
    end="Tmp/{sample}_select1ok.sam"
  input:
    "Tmp/{sample}.bam"
  log:
    "Logs/{sample}_extract.log"
  conda:
    "samtool.yaml"
  shell:
    "samtools view -H {input} > {output.select} ; "
    "cp {output.select} {output.end} "


rule hisat2_mapping:
  output:
    "Tmp/{sample}.bam"
  input:
#    R1=config["dataDir"]+"{sample}.fastq.gz" if config["rnaType"]=="single-end" else ??
    expand("Tmp/GenomeIdx/GenomeIdx.{idx}.ht2", idx=IDX),
    R1="Tmp/{sample}_R1.fastq",
    R2="Tmp/{sample}_R2.fastq"
  log:
    "Logs/{sample}_bwt2_mapping.log"
  conda: 
    "mapping.yml"
  shell:
    "hisat2 -x Tmp/GenomeIdx/GenomeIdx -k 1 --no-mixed --rna-strandness FR -1 {input.R1} -2 {input.R2} 2> {log} | samtools view -b -o {output}"
#  run: 
#    """
#     if (config["rnaType"]=="paired-end"):
#        shell("hisat2 -x Tmp/GenomeIdx -1 {input.R1} -2 {input.R2} -S {output} 2> {log}")
#     else:
#        shell("hisat2 -x Tmp/GenomeIdx -U {input.R1} -S {output} 2> {log}")
#    """
 

rule genome_hisat2_index:
  output:
    #directory("Tmp/GenomeIdx"),
    expand("Tmp/GenomeIdx/GenomeIdx.{ext}.ht2", ext=IDX)
  input:
    fna=config["dataDir"]+config["genome"]
  params:
    idx="Tmp/GenomeIdx/GenomeIdx"
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

