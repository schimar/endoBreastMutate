rule qc:
    input:
        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        r2 = lambda wildcards: getFqHome(wildcards.sample)[1]
    output:
        touch("raw/qc/fastqc/{sample}_qc.done")
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
        #r1 = "raw/{sample}_{lane}_R1_001.fastq.gz",
        #r2 = "raw/{sample}_{lane}_R2_001.fastq.gz",
        #r1 = expand("raw/{sample}_R1_001.fastq.gz", sample=sample_names),
        #r2 = expand("raw/{sample}_R2_001.fastq.gz", sample=sample_names),
        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
        adapters = config["adapters"]
    output:
        r1trmd = "trm/{sample}_R1.fq.gz",
        r2trmd = "trm/{sample}_R2.fq.gz"
    threads: 2
    message: """--- Quality trimming of fastq files."""
    shell: 
        """
        /home/schimar/bio/bbmap/bbduk.sh -Xmx1g in1={input.r1} in2={input.r2} out1={output.r1trmd} out2={output.r2trmd} trimq=6 qtrim=r hdist=1 bhist=trm/hist/{wildcards.sample}.bhist qhist=trm/hist/{wildcards.sample}.qhist lhist=trm/hist/{wildcards.sample}.lhist tpe tbo 
        """


rule clumpify:
    input:
        r1 = "trm/{sample}_{lane}_R1.fq.gz",
        r2 = "trm/{sample}_{lane}_R2.fq.gz"
    output:
        r1 = "trm/ddp/{sample}_{lane}_R1.fq.gz",
        r2 = "trm/ddp/{sample}_{lane}_R2.fq.gz",
    log: "trm/ddp/{sample}_{lane}.log"
    threads: 8
    message: """Deduplication of reads using clumpify."""
    shell:
        """
        /home/schimar/bio/bbmap/clumpify.sh Xmx16g in1={input.r1} in2={input.r2} out1={output.r1} out2={output.r2} dedupe=t ecc=t optical=t spany=t adjacent=t dupedist=40 unpair=t repair=t 2> {log}  
        """
# >> {output.smp} 2>&1

rule tadpole:
    input: 
        r1 = "trm/ddp/{sample}_{lane}_R1.fq.gz",
        r2 = "trm/ddp/{sample}_{lane}_R2.fq.gz"
    output:
        r1eco = "trm/eco/{sample}_{lane}_R1.fq.gz",
        r2eco = "trm/eco/{sample}_{lane}_R2.fq.gz",
    log: "trm/eco/{sample}_{lane}.log"
    priority: 50
    threads: 8
    message: """Error-correction of fastq files with tadpole."""
    shell:
        """
        /home/schimar/bio/bbmap/tadpole.sh -Xmx16g in1={input.r1} in2={input.r2} out1={output.r1eco} out2={output.r2eco} mode=correct k=50 overwrite=t 2> {log}
        """
# >> {output.smp} 2>&1

rule covStats:
    input:
        #
    output:
        # 
    threads: 8
    message: """Calculating coverage statistics."""

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
    message: """Concatenating fastq files for the same individuals."""
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





