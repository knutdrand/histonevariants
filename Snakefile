configfile:"config.json"
include: "repeat_analysis.sm"
include: "jung.sm"
track_hub = "../../var/www/html/trackhub_knut/mm10/"
NELS="u1452@nelstor0.cbu.uib.no:/elixir-chr/nels/users/u1452/Projects/UiO_Dahl_Chromatin_2018/MadeleineFosslie_MF/200110_A00943.B.Project_Fosslie-Libs12-2020-01-06/"
key="../u1452@nelstor0.cbu.uib.no.key"
repeat_names = "IAPEz-int  IAPLTR1a_Mm  IAPLTR1_Mm  L1Md_A  L1Md_F  L1Md_G L1Md_T  MTA_Mm  MTA_Mm-int  MTB".split()
my_names = config["chips"] #["1-p12-AZ-200", "2-p12-ac-100", "11-mor-AZ460", "12-mor-ac230"]
shorten_name = lambda name: ".".join(name.split("-")[1:3])


r_figure = """cat {input} | Rscript -e "library(tidyverse); png('{output}'); %s; dev.off()" """

rule all:
    input:
        expand(track_hub+"V2-2019-{name}.bw", name=config["samples"])

rule all_2015:
    input:
        expand(track_hub+"2015-{name}.bw", name=config["2015_samples"])

rule all_repeats:
    input:
        expand("repeats/{name}.bed.gz", name=["L1Md_T", "L1Md_A", "L1Md_F", "L1Md_G"])

rule all_repeats_cov:
    input:
        expand("domain_flank_repeat_enrichment/{name}.png", name=repeat_names)
rule all_flankrepeats:
    input:
        expand("flank_repeats/{repeat}/V2-2019-{name}.txt", repeat=repeat_names, name=["1-p12-AZ-200", "2-p12-ac-100", "11-mor-AZ460", "12-mor-ac230"])

rule clean_repeats:
    input:
        "repeats/{name}.bed.gz"
    output:
        "clean_repeats/{name}.bed.gz"
    shell:
        "zgrep -P 'chr[\d+,X, Y, M]' {input} | gzip > {output}"

rule all_fragment_sizes:
    input:
        expand("fragment_sizes/V2-2019-{name}.npy", name=config["samples"]),
        expand("fragment_sizes/V2-2019-{name}.png", name=config["samples"])

rule import_data:
    output:
        temp("reads/2019/{sample}_L{lane}_R{read}.fastq.gz")
    shell:
        "scp -i {key} {NELS}Sample_{wildcards.sample}/{wildcards.sample}_S*_L00{wildcards.lane}_R{wildcards.read}_001.fastq.gz {output}"

rule merge_lanes:
    input:
        "reads/2019/{sample}_L1_R{read}.fastq.gz",
        "reads/2019/{sample}_L2_R{read}.fastq.gz",
    output:
        temp("merged_reads/2019-{sample}_R{read}.fastq.gz")
    shell:
        "cat {input} > {output}"

rule merge_lanes_2015:
    input:
        "reads/2015/{sample}_L1_R{read}.fastq.gz",
        "reads/2015/{sample}_L2_R{read}.fastq.gz",
        "reads/2015/{sample}_L3_R{read}.fastq.gz",
        "reads/2015/{sample}_L4_R{read}.fastq.gz"
    output:
        temp("merged_reads/2015-{sample}_R{read}.fastq.gz")
    shell:
        "cat {input} > {output}"

rule trim_adaptors:
    input:
        "merged_reads/{name}_R1.fastq.gz",
        "merged_reads/{name}_R2.fastq.gz"
    output:
        temp("trimmed_reads/V2-{name}_R1.fastq.gz"),
        temp("trimmed_reads/V2-{name}_R2.fastq.gz"),
        'logs/cutadapt_{name}.log'
    threads: 16
    shell:
        'cutadapt --nextseq-trim=20 -m 10 -j {threads} -a "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC" -A "AATGATACGGCGACCACCGAGATCTACAC" -o {output[0]} -p {output[1]} {input} > {output[2]}'

rule bwa_map:
    input:
        "trimmed_reads/{name}_R1.fastq.gz",
        "trimmed_reads/{name}_R2.fastq.gz"
    output:
        temp("mapped_reads/{name}.bam")
    threads: 16
    shell:
        "bwa mem -t {threads} ../data/mm10.fa {input} | samtools view -Sb - > {output}"

rule filter_alignments:
    input:
        "mapped_reads/{name}.bam"
    output:
        "logs/{name}.flagstat",
        temp("filtered_alignments/{name}.bam")
    shell:
        """
	samtools flagstat {input} > {output[0]}
	samtools view -F 1804 -f 2 -u {input} > {output[1]}
	"""

rule filter_alignments_mapq:
    input:
        "filtered_alignments/{name}.bam"
    output:
        temp("filtered_alignments_mapq/{name}.bam")
    shell:
        "samtools view -q 30 -u {input} > {output}"

rule qc:
    input:
        "reads/{sample}.fastq.gz"
    output:
        "qc/{sample}_fastqc/fastqc_data.txt"
    shell:
        "fastqc {input} "

rule sort_bed:
    input:
        "{folder}/{sample}.bed"
    output:
        temp("{folder}_sorted/{sample}.bed")
    shell:
        "bedtools sort -i {input} > {output}"

rule fragment_bed:
    input:
        "filtered_alignments_mapq/{sample}.bam"
    output:
        temp("fragments/{sample}.bed")
    shell:
        "/usr/local/bin/macs2 randsample -i {input} -f BAMPE -p 100 -o {output}"

rule call_peaks:
    input:
        "repeat_fragments_small_removed_unique/{name}.bed.gz",
        lambda wildcards: "repeat_fragments_small_removed_unique/V2-2019-%s.bed.gz" % config["inputs"][wildcards.name.split("-")[3]]
    output:
        "macs_output/{name}_peaks.broadPeak",
        "macs_output/{name}_treat_pileup.bdg",
	"macs_output/{name}_control_lambda.bdg"
    shell:
        "macs2 callpeak -f BEDPE -t {input[0]} -c {input[1]} --bdg -n {wildcards.name} --broad --outdir macs_output/"
rule repeat_coverage_hack:
    input:
        "clean_repeats/{sample}.bed.gz"
    output:
        "repeat_coverage/{sample}.bdg"
    shell:
        "bedtools genomecov -bga -i {input} -g data/chrom.sizes.txt > {output}"

rule get_coverage:
    input:
        "fragments_unique/{sample}.bed"
    output:
        temp("coverage/{sample}.bdg")
    shell:
        "bedtools genomecov -bg -i {input} -g data/mm10.chrom.sizes > {output}"

rule filter_duplicates:
    input:
        "{folder}_sorted/{sample}.bed"
    output:
        "{folder}_unique/{sample}.bed"
    shell:
        "chiptools filterdup {input} > {output}"

rule make_track:
    input:
        "coverage/{sample}.bdg"
    output:
        temp("coverage/{sample}.bw")
    shell:
        "./bdg2bw {input} data/mm10.chrom.sizes"

rule get_fragment_sizes:
    input:
        "fragments_unique/{sample}.bed.gz"
    output:
        "fragment_sizes/{sample}.png",
        "fragment_sizes/{sample}.npy"
    shell:
        "zcat {input} | chiptools sizehist - 500 2 {output}"

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
	
rule filter_repeats:
    input:
        "data/mouse_repeats.bed.gz"
    output:
        "repeats/{name}.bed.gz"
    shell:
        "zcat {input} | awk '{{if (($3-$2)>250){{print}}}}' | grep {wildcards.name} | gzip > {output}"

rule import_antibody_data:
    output:
        "reads/2015/{sample}_{type}_L{lane}_R{read}.fastq.gz"
    shell:
        "mv 2015_reads/{wildcards.sample}_Z1-{wildcards.type}_L00{wildcards.lane}_R{wildcards.read}_001.fastq.gz {output}"

rule get_directed_flanks:
    input:
        "2016/2016_mouse_Oocyte_H3K4me3_domains_merged.bed",
    output:
        temp("2016/M_borders+.bed"),
        temp("2016/M_borders-.bed"),
        "2016/M_borders.bed"
    shell:
        """
        awk '{{OFS="\t"}}{{print $1, $2, $3, ".", ".", "+"}}' {input} > {output[0]}
        awk '{{OFS="\t"}}{{print $1, $2, $3, ".", ".", "-"}}' {input} >> {output[1]}
        paste -d '\n' {output[0]} {output[1]} > {output[2]}
        """
#         """
#         bedtools flank -i {input[0]} -g {input[1]} -l 1000 -r 0| awk '{{OFS="\t"}}{{print $1, $2, $3, ".", ".", "+"}}' > {output}
#         bedtools flank -i {input[0]} -g {input[1]} -r 1000 -l 0| awk '{{OFS="\t"}}{{print $1, $2, $3, ".", ".", "-"}}' >> {output}
#         """

rule get_flank_enrichment:
    input:
        "regions/non_tss_containing_domains/M_borders.bed",
        "macs_output/{name}_treat_pileup.bdg"
    output:
        "domain_flank_enrichment/{name}.npy",
        "domain_flank_enrichment/{name}.png"
    shell:
        "cat {input[1]} | chiptools tssplot {input[0]} {output}"

rule get_flank_repeat_enrichment:
    input:
        "regions/non_tss_containing_domains/M_borders.bed",
        "repeat_coverage/{name}.bdg"
    output:
        "domain_flank_repeat_enrichment/{name}.npy",
        "domain_flank_repeat_enrichment/{name}.png"
    shell:
        "cat {input[1]} | chiptools tssplot {input[0]} {output}"        

rule download_chrom_sizes:
    output:
        temp("data/chromInfo.txt.gz")
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/chromInfo.txt.gz -O {output}"

rule clean_chrom_sizes:
    input:
        "data/chromInfo.txt.gz"
    output:
        "data/chrom.sizes.txt"
    shell:
        "zgrep -P 'chr[\d+,X, Y, M]' {input} > {output}"


rule get_tss:
    input:
        "data/genes.bed"
    output:
        "regions/tss/1.bed"
    shell:
        """awk '{{OFS="\t"}}{{if ($6=="+") {{print $1,$2,$2+1}} else {{print $1, $3-1, $3}}}}' {input} | uniq > {output}"""

rule get_tss_containing_domains:
    input:
        "2016/M_borders.bed",
        "regions/tss/1.bed"
    output:
        "regions/tss_containing_domains/M_borders.bed"
    shell:
        """bedtools intersect -a {input[0]} -b {input[1]} -c | awk '{{if ($NF>0) print}}' > {output}"""

rule get_non_tss_containing_domains:
    input:
        "2016/M_borders.bed",
        "regions/tss/1.bed"
    output:
        "regions/non_tss_containing_domains/M_borders.bed"
    shell:
        """bedtools intersect -a {input[0]} -b {input[1]} -c | awk '{{if ($NF==0) print}}' > {output}"""

rule get_flanks:
    input:
        "{name}.bed",
        "data/chrom.sizes.txt"
    output:
        "{name}.flanks{b}.bed"
    shell:
        """bedtools flank -i {input[0]} -b {wildcards.b} -g {input[1]} > {output}"""

rule flank_repeat_enrichment_analysis:
    input:
        "2016/2016_mouse_Oocyte_H3K4me3_domains_merged.flank.bed",
        "clean_repeats/{repeat}.bed.gz",
        "macs_output/{name}_treat_pileup.bdg.gz",
        "data/chrom.sizes.txt"
    output:
        "flank_repeats/{repeat}/{name}.txt"
    shell:
        """
        sortBed -i {input[1]} -g {input[3]} | intersectBed -a stdin -b {input[0]} | intersectBed -a {input[2]} -b stdin | awk '{{v+=($3-$2)*$4;n+=($3-$2)}}END{{print v, n, (v/n)}}'> {output}
        sortBed -i {input[1]} -g {input[3]} | complementBed -i stdin -g {input[3]} | intersectBed -a stdin -b {input[0]} | intersectBed -a {input[2]} -b stdin | awk '{{v+=($3-$2)*$4;n+=($3-$2)}}END{{print v, n, (v/n)}}'>> {output}
        """

rule calculate_enrichment:
    input:
        "flank_repeats/{repeat}/{name}.txt"
    output:
        "flank_repeats/{repeat}/{name}.enrich"
    shell:
        """cat {input} | python -c "import sys; tmp=[float(line.split()[-1]) for line in sys.stdin];print(tmp[0]/tmp[1])" > {output}"""               

rule summarize_folder:
    input:
        expand("flank_repeats/{{repeat}}/V2-2019-{name}.enrich", name=my_names)
    output:
        "flank_repeats/{repeat}/ALL.csv"
    shell:
        "paste -d ',' {input} > {output}"

rule summarize_all_repeats:
    input:
        expand("flank_repeats/{repeat}/ALL.csv", repeat=repeat_names)
    output:
        temp("flank_repeats/header.csv"),
        temp("flank_repeats/data.csv"),
        "flank_repeats/ALL.csv"
    shell:
        """
        echo "%s" > {output[0]}
        echo "%s" > {output[1]}
        cat {input} >> {output[1]}
        paste -d ',' {output[0]} {output[1]} > {output[2]}
        """%("\n".join(["Repeat"] +repeat_names), ",".join(shorten_name(n) for n in my_names))

rule get_bdg_average:
    input:
        "{name}.bdg.gz"
    output:
        "{name}.bdg.average"
    shell:
        "zcat {input} | awk '{{v+=($3-$2)*$4;n+=($3-$2)}}END{{print v, n, (v/n)}}' > {output}"

rule heatmap:
    input:
        "{name}.csv"
    output:
        "{name}.heatmap.png"
    shell:
        r_figure % "read.csv('stdin', row.names=1) %>% as.matrix %>% heatmap(scale='none')"
