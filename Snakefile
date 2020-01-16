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
        "trimmed_reads/{name}_R1_001.fastq.gz",
        "trimmed_reads/{name}_R2_001.fastq.gz"
    output:
        "mapped_reads/{name}.bam"
    threads: 16
    shell:
        "bwa mem -t {threads} ../data/mm10.fa {input} | samtools view -Sb - > {output}"

rule filter_alignments:
    input:
        "mapped_reads/{name}.bam"
    output:
        "logs/{name}.flagstat",
        "filtered_alignments/{name}.bam"
    shell:
        """
	samtools flagstat {input} > {output[0]}
	samtools view -F 1804 -f 2 -u {input} > {output[1]}
	"""

rule qc:
    input:
        "reads/{sample}.fastq.gz"
    output:
        "qc/{sample}_fastqc/fastqc_data.txt"
    shell:
        "fastqc {input}"

rule sort_bed:
    input:
        "fragments/{sample}.bed"
    output:
        "sorted_fragments/{sample}.bed"
    shell:
        "bedtools sort -i {input} > {output}"

rule fragment_bed:
    input:
        "filtered_alignments/{sample}.bam"
    output:
        temp("fragments/{sample}.bed")
    shell:
        "/usr/local/bin/macs2 randsample -i {input} -f BAMPE -p 100 -o {output}"

rule get_coverage:
    input:
        "unique_fragments/{sample}.bed"
    output:
        temp("coverage/{sample}.bdg")
    shell:
        "bedtools genomecov -bg -i {input} -g data/mm10.chrom.sizes > {output}"

rule filter_duplicates:
    input:
        "sorted_fragments/{sample}.bed"
    output:
        "unique_fragments/{sample}.bed"
    shell:
        "macs2 filterdup -i {input} --keep-dup=1 -o {output}"

rule make_track:
    input:
        "coverage/{sample}.bdg"
    output:
        temp("coverage/{sample}.bw")
    shell:
        "./bdg2bw {input} data/mm10.chrom.sizes"

rule export_track:
    input:
        "coverage/{sample}.bw"
    output:
        "/var/www/html/trackhub_knut/mm10/{sample}.bw"
    shell:
        "mv {input} {output}"
