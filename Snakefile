rule trim_adaptors:
    input:
        "reads/{name}_R1_001.fastq.gz",
        "reads/{name}_R2_001.fastq.gz"
    output:
        temp("trimmed_reads/{name}_R1_001.fastq.gz"),
        temp("trimmed_reads/{name}_R2_001.fastq.gz")
    shell:
        'cutadapt -a "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" -A "AATGATACGGCGACCACCGAGATCTACAC" -o {output[0]} -p {output[1]} {input}'

rule bwa_map:
    input:
        "../data/mm10.fa",
        "trimmed_reads/{name}_R1_001.fastq.gz",
        "trimmed_reads/{name}_R2_001.fastq.gz"
    output:
        "mapped_reads/{name}.bam"
    threads: 16
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

rule qc:
    input:
        "reads/{sample}.fastq.gz"
    output:
        "qc/{sample}_fastqc/fastqc_data.txt"
    shell:
        "fastqc {input}"

rule fragment_bed:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "fragment_files/{sample}_fragments.bed"
    shell:
        "macs2 randsample -i {input} -f BAMPE -p 100 -o {output}"
