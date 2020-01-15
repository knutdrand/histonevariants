rule bwa_map:
    input:
        "data/genome.fa",
	"reads/{name}_R1_001.fastq.gz",
        "reads/{name}_R2_001.fastq.gz"
    output:
        "mapped_reads/A.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

rule qc:
    input:
        "reads/{sample}.fastq.gz"
    output:
        "qc/{sample}_fastqc/fastqc_data.txt"
    shell:
        "fastqc {input}"
