SAMPLES, = glob_wildcards(config["dataDir"]+"{sample}_R1.fastq.gz")
IDX = ["1","2","3","4","5","6","7","8"]

rule all :
  input :
    expand("FastQC/{sample}_R1_fastqc.html", sample=SAMPLES),
    expand("FastQC/{sample}_R2_fastqc.html", sample=SAMPLES),
    expand(config["dataRes"]+"{sample}_"+config["chr"]+"_"+str(config["nbAlignRandom"])+"r_"+str(config["nbAlignChr"])+"c_R{R}.fastq.gz", sample=SAMPLES, R=["1","2"]),


rule compress:
  output:
    config["dataRes"]+"{sampleRed}"+config["chr"]+"_"+str(config["nbAlignRandom"])+"r_"+str(config["nbAlignChr"])+"c_R{R}.fastq.gz",
  input:
    config["dataRes"]+"{sampleRed}"+config["chr"]+"_"+str(config["nbAlignRandom"])+"r_"+str(config["nbAlignChr"])+"c_R{R}.fastq",
  log:
    out="Logs/compress_{sampleRed}"+config["chr"]+"_"+str(config["nbAlignRandom"])+"r_"+str(config["nbAlignChr"])+"c_R{R}.out",
    err="Logs/compress_{sampleRed}"+config["chr"]+"_"+str(config["nbAlignRandom"])+"r_"+str(config["nbAlignChr"])+"c_R{R}.err",
  shell:
    """
    gzip {input} 2>{log.err}
    """


rule create_fastq:
  output:
    R1=config["dataRes"]+"{sample}_"+config["chr"]+"_"+str(config["nbAlignRandom"])+"r_"+str(config["nbAlignChr"])+"c_R1.fastq",
    R2=config["dataRes"]+"{sample}_"+config["chr"]+"_"+str(config["nbAlignRandom"])+"r_"+str(config["nbAlignChr"])+"c_R2.fastq",
  input:
    "Tmp/{sample}_select.sam",
  log:
    out="Logs/{sample}_createFastqRed.log",
    err="Logs/{sample}_createFastqRed.err",
  conda: config["condaDir"]+"ce_RNASeqReduction.yml",
  envmodules: "samtools",
  shell:
    """
    samtools fastq -n -1 {output.R1} -2 {output.R2} Tmp/{wildcards.sample}_select.sam 1>{log.out} 2>{log.err} 
    """

    
rule extract_reads:
  output:
    temp("Tmp/{sample}_select.sam"),
  input:
    "Tmp/{sample}.bam",
  log:
    out="Logs/{sample}_extract.our",
    err="Logs/{sample}_extract.err",
  conda: config["condaDir"]+"ce_RNASeqReduction.yml",
  envmodules: "samtools",
  shell:
    "samtools view -H {input} > {output} 2>{log.err} ; "
    "set +o pipefail ; "
    "samtools view {input} | head -n "+str(config["nbAlignRandom"])+" | grep -v "+config["chr"]+" >> {output} 2>>{log.err} ; "
    "samtools view {input} | grep "+config["chr"]+" | head -n "+str(config["nbAlignChr"])+" >> {output} 2>>{log.err} ; "


rule hisat2_mapping:
  output:
    temp("Tmp/{sample}.bam"),
  input:
    expand("Tmp/GenomeIdx/GenomeIdx.{idx}.ht2", idx=IDX),
    R1="Tmp/{sample}_R1.fastq",
    R2="Tmp/{sample}_R2.fastq",
  params:
    sepe=config["rnaType"],
  threads: config["threads_mapping"]
  resources:
    mem_mb=config["mem_mapping"]
  log:
    out="Logs/{sample}_bwt2_mapping.out",
    err="Logs/{sample}_bwt2_mapping.err",
  conda: config["condaDir"]+"ce_RNASeqReduction.yml",
  envmodules: 
     "hisat2", 
     "samtools",
  shell:
    """
    if [ {params.sepe} = "paire-end" ]; then 
       hisat2 -x Tmp/GenomeIdx/GenomeIdx -p {threads} -k 1 --no-mixed --rna-strandness FR -1 {input.R1} -2 {input.R2} 2>{log.err} | samtools view -@ {threads} -b -o {output} 2>>{log.err}
    else # single-end case:
       hisat2 -x Tmp/GenomeIdx/GenomeIdx -p {threads} -k 1 -U {input.R1} 2>{log.err} | samtools view -@ {threads} -b -o {output} 2>>{log.err}
    fi
    """
 

rule genome_hisat2_index:
  output:
    temp(expand("Tmp/GenomeIdx/GenomeIdx.{ext}.ht2", ext=IDX)),
  input:
    fna=config["genome"],
  params:
    idx="Tmp/GenomeIdx/GenomeIdx",
  resources:
    mem_mb=config["mem_index"],
  conda: config["condaDir"]+"ce_RNASeqReduction.yml",
  envmodules: "hisat2",
  log:
    out="Logs/genome_hisat2_index.out",
    err="Logs/genome_hisat2_index.err",
  shell: 
    "hisat2-build {input.fna} {params.idx} 1>{log.out} 2>{log.err}"


rule fastqc:
  output: 
    zip="FastQC/{sample}_fastqc.zip",
    html="FastQC/{sample}_fastqc.html",
  input: 
    fastq = "Tmp/{sample}.fastq",
    previous_gunzip_ok = "Tmp/{sample}.uncompress.done",
  params:
    sepe=config["rnaType"],
  log:
    out="Logs/{sample}_fastqc.out",
    err="Logs/{sample}_fastqc.err",
  conda: config["condaDir"]+"ce_RNASeqReduction.yml",
  envmodules: "fastqc",
  shell: 
    """
    if [ {params.sepe} = "single-end" ]; then 
       touch {output.zip} {output.html}
    else    
       fastqc --outdir FastQC/ {input.fastq} 1>{log.out} 2>{log.err}
    fi
    """


rule uncompress:
  output:
    fastq = temp("Tmp/{sample}.fastq"),
    gunzip_done = temp(touch("Tmp/{sample}.uncompress.done")), # end check to avoid latency error
  input:
    config["dataDir"]+"{sample}.fastq.gz",
  params:
    sepe=config["rnaType"],
  log:
    err="Logs/{sample}_uncompress.err",
  shell: 
    """
    gunzip -c {input} > {output.fastq} 2>{log.err} ;
    # single-end case:
    if [ {params.sepe} = "single-end" ]; then 
       set +o pipefail ;
       touch `echo {output.fastq} | awk '{{gsub("R1.fastq","R2.fastq",$0);print $0}}' ` 
    fi  
    """


