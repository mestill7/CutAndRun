SAMPLES, = glob_wildcards("fastq/{sample}_R1_001.fastq.gz")
configfile: "config.yaml"

ALL_BW = expand('processed/bw/{sample}.unique.sorted.rmdup.chr.bw',
                        sample=SAMPLES)
ALL_BED = expand('processed/bed/{sample}.unique.sorted.rmdup.chr.bed',
                        sample=SAMPLES)
ALL_BAM = expand('processed/bam/{sample}.unique.sorted.rmdup.bam',
                        sample=SAMPLES)

COUNTS_MATRIX = "processed/htseq_counts_matrix.txt"
MULTIQC_REPORT = "multiqc_report.html"

if config["experiment"] == "cutrun":
    rule all:
        input: ALL_BED, ALL_BW, ALL_BAM, MULTIQC_REPORT

    rule trim_fastq_fastqc:
        input:
            pair1 = "fastq/{sample}_R1_001.fastq.gz",
            pair2 = "fastq/{sample}_R2_001.fastq.gz"
        output:
            trimmed_pair1 = temp("logs/{sample}_R1_001_val_1.fq.gz"),
            trimmed_pair2 = temp("logs/{sample}_R2_001_val_2.fq.gz"),
            fastqc_zipfile1 = "fastqc/{sample}_R1_001_fastqc.zip",
            fastqc_zipfile2 = "fastqc/{sample}_R2_001_fastqc.zip"
        log:
            "logs/{sample}.trim_adapters.log"
        run:
            shell("trim_galore {input.pair1} {input.pair2} --paired -o ./logs")
            shell("fastqc {input.pair1} {input.pair2} -o ./fastqc")

    rule split_length_long:
        input:
            trimg_pair1 = "logs/{sample}_R1_001_val_1.fq.gz",
            trimg_pair2 = "logs/{sample}_R2_001_val_2.fq.gz"
        output:
            cut_r1_p1 = temp("logs/{sample}_t1_R1.len75.fastq"),
            cut_r2_p1 = temp("logs/{sample}_t1_R2.len75.fastq")
        params:
            read_length = config["read_length"]
        log:
            "logs/{sample}.split_length_keeplong.log"
        run:
            shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r1_p1} {input.trimg_pair1}"),
            shell("cutadapt --minimum-length {params.read_length} -o {output.cut_r2_p1} {input.trimg_pair2}")

    rule split_length_short:
        input:
            trimg_pair1 = "logs/{sample}_R1_001_val_1.fq.gz",
            trimg_pair2 = "logs/{sample}_R2_001_val_2.fq.gz"
        output:
            cut_r1_p2 = temp("logs/{sample}_t1_R1.lt75.fastq"),
            cut_r2_p2 = temp("logs/{sample}_t1_R2.lt75.fastq")
        params:
            read_length = config["read_length_max"]
        log:
            "logs/{sample}.split_length_keepshort.log"
        run:
            shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r1_p2} {input.trimg_pair1}"),
            shell("cutadapt --maximum-length {params.read_length} -o {output.cut_r2_p2} {input.trimg_pair2}")

    rule trim_long:
        input:
            cut_r1_p1 = "logs/{sample}_t1_R1.len75.fastq",
            cut_r2_p1 = "logs/{sample}_t1_R2.len75.fastq"
        output:
            cut_r1_p3 = temp("logs/{sample}_t1_R1.len75_trim.fastq"),
            cut_r2_p3 = temp("logs/{sample}_t1_R2.len75_trim.fastq")
        log:
            "logs/{sample}.split_length_keepshort.log"
        run:
            shell("cutadapt -u -6 -o {output.cut_r1_p3} {input.cut_r1_p1}"),
            shell("cutadapt -u -6 -o {output.cut_r2_p3} {input.cut_r2_p1}")

    rule combine_split_lengths:
        input:
            cut_r1_p3 = "logs/{sample}_t1_R1.len75_trim.fastq",
            cut_r2_p3 = "logs/{sample}_t1_R2.len75_trim.fastq",
            cut_r1_p2 = "logs/{sample}_t1_R1.lt75.fastq",
            cut_r2_p2 = "logs/{sample}_t1_R2.lt75.fastq"
        output:
            cut_r1_p4 = temp("logs/{sample}_t2_R1.fastq"),
            cut_r2_p4 = temp("logs/{sample}_t2_R2.fastq")
        run:
            shell("cat {input.cut_r1_p3} {input.cut_r1_p2} > {output.cut_r1_p4}"),
            shell("cat {input.cut_r2_p3} {input.cut_r2_p2} > {output.cut_r2_p4}")

    rule sort_combined_lengths:
        input:
            cut_r1_p4 = "logs/{sample}_t2_R1.fastq",
            cut_r2_p4 = "logs/{sample}_t2_R2.fastq"
        output:
            cut_r1_p5 = temp("logs/{sample}_t2_R1_sorted.fastq.gz"),
            cut_r2_p5 = temp("logs/{sample}_t2_R2_sorted.fastq.gz")
        run:
            shell("module load java"),
            shell("/sc/hydra/work/estilm01/Tools/ngsutilsj-master/dist/ngsutilsj fastq-sort --output logs/{wildcards.sample}_t2_R1_sorted.fastq.gz {input.cut_r1_p4}"),
            shell("/sc/hydra/work/estilm01/Tools/ngsutilsj-master/dist/ngsutilsj fastq-sort --output logs/{wildcards.sample}_t2_R2_sorted.fastq.gz {input.cut_r2_p4}")

    rule bowtie2:
        input:
            trimmed_pair1 = "logs/{sample}_t2_R2_sorted.fastq.gz",
            trimmed_pair2 = "logs/{sample}_t2_R2_sorted.fastq.gz"
        params:
            index = config["index"]
        output:
            bam = temp("processed/bam/{sample}.sorted.bam"),
            bambai = temp("processed/bam/{sample}.sorted.bam.bai")
        threads:
            config["threads_for_alignment"]
        log:
            "logs/{sample}.alignment.log"
        run:
            shell("bowtie2 -p {threads} --dovetail --phred33 -x {params.index} -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} 2> {log} > processed/bam/{wildcards.sample}.sam"),
            shell("samtools sort processed/bam/{wildcards.sample}.sam | samtools view -bS - > processed/bam/{wildcards.sample}.bam"),
            shell("rm processed/bam/{wildcards.sample}.sam"),
            shell("samtools view -bh -f 3 -F 4 -F 8 processed/bam/{wildcards.sample}.bam > processed/bam/{wildcards.sample}_mapped.bam"),
            shell("samtools index processed/bam/{wildcards.sample}_mapped.bam"),
            shell("samtools sort processed/bam/{wildcards.sample}_mapped.bam > {output.bam}"),
            shell("samtools index {output.bam}")
            shell("rm processed/bam/{wildcards.sample}_mapped.bam*")

    rule rmdup:
        input:
            "processed/bam/{sample}.sorted.bam"
        output:
            bam = "processed/bam/{sample}.unique.sorted.rmdup.bam"
        params:
            picardbin = "/hpc/packages/minerva-common/picard/2.7.1",
            picardmetric = "logs/{sample}.markdups.metrics.txt"
        run:
            shell("picard MarkDuplicates INPUT={input} OUTPUT={output.bam} VALIDATION_STRINGENCY=SILENT METRICS_FILE={params.picardmetric}"),
            shell("samtools index {output.bam}")

    rule rmdup_to_chrbam:
        input:
            dup_removed = "processed/bam/{sample}.unique.sorted.rmdup.bam"
        output:
            chrbam = "processed/bam/{sample}.unique.sorted.rmdup.chr.bam"
        log:
            "logs/{sample}.chrbam.log"
        shell:
            'samtools view -H {input.dup_removed} '
            '| sed -e "s/SN:\([0-9XY]\)/SN:chr\\1/" -e "s/SN:MT/SN:chrM/" '
            '| samtools reheader - {input.dup_removed} > {output.chrbam}'

    rule chrbam_to_bw:
        input:
            chrbam = "processed/bam/{sample}.unique.sorted.rmdup.chr.bam"
        output:
            bw_file = "processed/bw/{sample}.unique.sorted.rmdup.chr.bw"
        log:
            "logs/{sample}.bw.log"
        run:
            shell("samtools index {input.chrbam}"),
            shell("bamCoverage -b {input.chrbam} -o {output.bw_file} --binSize 10 --normalizeUsingRPKM")

    rule chrbam_to_bed:
        input:
            chrbam = "processed/bam/{sample}.unique.sorted.rmdup.chr.bam"
        output:
            bed = "processed/bed/{sample}.unique.sorted.rmdup.chr.bed"
        log:
            "logs/{sample}.bed.log"
        shell:
            "bedtools bamtobed -i {input.chrbam} > {output.bed} 2> {log}"

    rule run_multiqc:
        input:
            bed = expand('processed/bed/{sample}.unique.sorted.rmdup.chr.bed',sample=SAMPLES)
        output:
            multiqc_report = "multiqc_report.html"
        params:
            multiqc_config = config["multiqc_yaml"]
        shell:
            "multiqc . -f --config {params.multiqc_config}"
