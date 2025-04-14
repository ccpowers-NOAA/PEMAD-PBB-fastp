import os

samples = []
for file in os.listdir("00_reads"):
    if file.endswith("R1.fastq.gz"):
        samples.append(file)

rule all:
    input:
        "log/multiqc.done"

rule pre_fastqc:
    output:
        "log/pre_fastqc_{sample}.done"
    conda:
        "envs/fastp-0.24.0.yaml"
    resources:
        runtime=1000,
        cpus_per_task=1,
        mem_mb=1000
    shell:
        """
        mkdir -p log
        mkdir -p pretrimming
        fastqc 00_reads/{wildcards.sample}
        filename={wildcards.sample}
        fastqc 00_reads/${{filename/R1.fastq.gz/R2.fastq.gz}}
        mv 00_reads/*html pretrimming
        mv 00_reads/*html pretrimming
        touch {output}
        """

rule fastp:
    input:
        "log/pre_fastqc_{sample}.done"
    output:
        "01_trimmed/{sample}_trimmed.fastq.gz"
    conda:
        "envs/fastp-0.24.0.yaml"
    resources:
        runtime=1000,
        cpus_per_task=1,
        mem_mb=1000
    shell:
        """
        mkdir -p 01_trimmed
        filename={wildcards.sample}
        fastp -i 00_reads/{wildcards.sample} -I 00_reads/${{filename/R1.fastq.gz/R2.fastq.gz}} \
          -o 01_trimmed/{wildcards.sample}_trimmed.fastq.gz -O 01_trimmed/${{filename/R1.fastq.gz/R2.fastq.gz}}_trimmed.fastq.gz
        """


rule post_fastqc:
    input:
        "01_trimmed/{sample}_trimmed.fastq.gz"
    output:
        "log/post_fastqc_{sample}.done"
    conda:
        "envs/fastp-0.24.0.yaml"
    resources:
        runtime=1000,
        cpus_per_task=1,
        mem_mb=1000
    shell:
        """
        mkdir -p posttrimming
        fastqc 01_trimmed/{wildcards.sample}_trimmed.fastq.gz
        filename={wildcards.sample}
        fastqc 01_trimmed/${{filename/R1.fastq.gz/R2.fastq.gz}}_trimmed.fastq.gz
        mv 00_reads/*html posttrimming
        mv 00_reads/*html posttrimming
        touch {output}
        """

rule multiqc:
    input:
        expand("log/post_fastqc_{sample}.done", sample=samples)
    output:
        "log/multiqc.done"
    conda:
        "envs/fastp-0.24.0.yaml"
    resources:
        runtime=1000,
        cpus_per_task=1,
        mem_mb=1000
    shell:
        """
        multiqc pretrimming -n pretrimming.html
        multiqc posttrimming -n posttrimming.html
        touch {output}
        """
