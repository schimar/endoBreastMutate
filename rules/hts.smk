rule qc:
    input:
        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        r2 = lambda wildcards: getFqHome(wildcards.sample)[1]
    output:
        touch("raw/qc/fastqc/{sample}_qc.done")
        #"raw/qc/{sample}_R1_001_fastqc.html",
        #"raw/qc/{sample}_R1_001_fastqc.zip",
        #"raw/qc/{sample}_R2_001_fastqc.html",
        #"raw/qc/{sample}_R2_001_fastqc.zip"
    #log: "log/{sample}.qc.log.txt"
    #resources:
    #    mem = 1000,
    #    time = 300
    threads: 1
    message: """--- Quality check of raw data with FastQC before trimming."""
    shell: """
        fastqc -o raw/qc/fastqc/ -f fastq {input.r1} &
        fastqc -o raw/qc/fastqc/ -f fastq {input.r2}
        """


rule trim:
    input:
        #r1 = expand("raw/{sample}_R1_001.fastq.gz", sample=sample_names),
        #r2 = expand("raw/{sample}_R2_001.fastq.gz", sample=sample_names),
        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
        adapters = config["adapters"]
    output:
        r1trmd = "trm/{sample}_R1.fq.gz",
        r2trmd = "trm/{sample}_R2.fq.gz"
    threads: 2
    message: """--- Quality trimming of fastq files before mapping."""
    shell: 
        """
        /home/schimar/bio/bbmap/bbduk.sh -Xmx1g in1={input.r1} in2={input.r2} out1={output.r1trmd} out2={output.r2trmd} trimq=6 qtrim=r hdist=1 bhist=trm/hist/{wildcards.sample}.bhist qhist=trm/hist/{wildcards.sample}.qhist lhist=trm/hist/{wildcards.sample}.lhist tpe tbo 
        """


rule mergeFQs:
    input:
        l1r1 = "trm/{idskeys}_L001_R1.fq.gz", 
        l1r2 = "trm/{idskeys}_L001_R2.fq.gz", 
        l2r1 = "trm/{idskeys}_L002_R1.fq.gz", 
        l2r2 = "trm/{idskeys}_L002_R2.fq.gz", 
        l3r1 = "trm/{idskeys}_L003_R1.fq.gz", 
        l3r2 = "trm/{idskeys}_L003_R2.fq.gz", 
        l4r1 = "trm/{idskeys}_L004_R1.fq.gz", 
        l4r2 = "trm/{idskeys}_L004_R2.fq.gz" 
    output:
        r1 = "trm/mrg/{idskeys}_R1.fq.gz",
        r2 = "trm/mrg/{idskeys}_R2.fq.gz" 
    shell:
        """
        cat {input.l1r1} {input.l2r1} {input.l3r1} {input.l4r1} > {output.r1} 
        cat {input.l1r2} {input.l2r2} {input.l3r2} {input.l4r2} > {output.r2} 
        """


      #  commands = [
      #      "cat trm/\* > trm/mrg/{output.R1}"    
      #      ] 
      #  for c in commands:
      #    shell(c)

 
rule refIndex:
	input:
		ref = config['ref']
	output:
		'ref/genome/1/summary.txt'
	shell:
		"""
		/home/schimar/bio/bbmap/bbmap.sh -Xmx24g ref={input.ref}
		"""


rule map:
	input:
		tr1 = lambda wildcards: getTrmHome(wildcards.sample)[0],
		tr2 = lambda wildcards: getTrmHome(wildcards.sample)[1],
		ref = config["ref"]
	output:
		bam = "map/{sample}.bam"
	threads: 12
	message: """--- Mapping reads to reference genome ---"""
	shell:
		"""
		/home/schimar/bio/bbmap/bbmap.sh -Xmx20g t={threads} ref={input.ref} in1={input.tr1} in2={input.tr2} out=stdout.sam minid=0.85 rgid={wildcards.sample} rglb=igaDNA rgsm={wildcards.sample} rgpl=ILLUMINA overwrite=f unpigz=t | samtools view -F 4 -Shu - | samtools sort - -o {output.bam}
		"""

rule bamIndex:
	input: 
		aln = 'map/{sample}.bam'	
	output:
		touch("map/{sample}.bamIndex.done")
	threads: 2
	message: """--- Indexing with samtools ---"""
	shell:
		"""
		samtools index {input.aln} 
		"""





