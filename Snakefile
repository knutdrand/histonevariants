configfile:"config.json"

track_hub = "../../var/www/html/trackhub_knut/mm10/"
NELS="u1452@nelstor0.cbu.uib.no:/elixir-chr/nels/users/u1452/Projects/UiO_Dahl_Chromatin_2018/MadeleineFosslie_MF/200110_A00943.B.Project_Fosslie-Libs12-2020-01-06/"
key="../u1452@nelstor0.cbu.uib.no.key"

rule all:
    input:
        expand(track_hub+"{name}.bw", name=config["samples"])

rule import_data:
    output:
        "reads/{sample}_L{lane}_R{read}.fastq.gz"
    shell:
        "scp -i {key} {NELS}Sample_{wildcards.sample}/{wildcards.sample}_S*_L00{wildcards.lane}_R{wildcards.read}_001.fastq.gz {output}"

rule merge_lanes:
    input:
        "reads/{sample}_L1_R{read}.fastq.gz",
        "reads/{sample}_L2_R{read}.fastq.gz"
    output:
        temp("merged_reads/{sample}_R{read}.fastq.gz")
    shell:
        "cat {input} > {output}"

rule trim_adaptors:
    input:
        "merged_reads/{name}_R1.fastq.gz",
        "merged_reads/{name}_R2.fastq.gz"
    output:
        temp("trimmed_reads/{name}_R1.fastq.gz"),
        temp("trimmed_reads/{name}_R2.fastq.gz")
    shell:
        'cutadapt -a "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" -A "AATGATACGGCGACCACCGAGATCTACAC" -o {output[0]} -p {output[1]} {input}'

rule bwa_map:
    input:
        "trimmed_reads/{name}_R1.fastq.gz",
        "trimmed_reads/{name}_R2.fastq.gz"
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
        track_hub+ "{sample}.bw"
    shell:
        """
	mv {input} {output}
	python3 ../create_hub.py ../../var/www/html/trackhub_knut/mm10/ > ../../var/www/html/trackhub_knut/mm10/trackDb.txt
	"""
	
