"""
This is a comment!

"""



rule fastp_pe:
    input:
        sample=["resources/downloads/{sample}_R1_001.fastq.gz", "resources/downloads/{sample}_R2_001.fastq.gz"]
    output:
        trimmed=["results/1.trimmed/{sample}_l.fastq", "results/1.trimmed/{sample}_r.fastq"],
        html="results/report/{sample}.html",
        json="results/report/{sample}.json"
    log:
        "results/logs/fastp/{sample}.log"
    params:
        extra=" "
    threads: 2
    shell:
        """
        fastp --in1 {input.sample[0]} \
              --in2 {input.sample[1]}  \
              --out1 {output.trimmed[0]} \
              --out2 {output.trimmed[1]} \
              --html {output.html} \
              --json {output.json}

        """
              




rule bowtie2:
    input:
        sample=["results/1.trimmed/{sample}_l.fastq", "results/1.trimmed/{sample}_r.fastq"]
    output:
        "results/2.mapped/{sample}.bam"
    log:
        "results/logs/bowtie2/{sample}.log"
    params:
        index="resources/refs/unaligned_refs_site_design.fasta",  # prefix of reference genome index (built with bowtie2-build)
        extra="--local --no-unal"  # optional parameters
    threads: 1
    shell:
        """
        bowtie2    {params.extra} \
                -x {params.index} \
                -1 {input.sample[0]} \
                -2 {input.sample[1]} \
                 | samtools view -hb - > {output}
      """
                



rule samtools_sort:
    input:
        "results/2.mapped/{sample}.bam"
    output:
        bam="results/2.mapped/{sample}.sorted.bam",
        bai="results/2.mapped/{sample}.sorted.bam.bai"
    shell:
        "samtools sort  {input} >{output.bam} ;"
        "samtools index {output.bam} " 



rule get_uniq_seqs:
    input:
        "results/2.mapped/{sample}.sorted.bam"
    params:
       script = 'workflow/scripts/uniq_reads.py'
    output:
        fasta="results/3.uniq_reads/{sample}_uniq.fasta"
    shell:
       """
       echo $( which python )
        $( which python ) {params.script}  -b {input} -o {output}

       """

rule find_mutations:
    input:
        "results/3.uniq_reads/{sample}_uniq.fasta"
    params:
        script = 'workflow/scripts/batch_mutation_search.py'
    output:
        "results/4.tsr_report/{sample}_tsr_report.tsv"
    shell:
       """
        $( which python ) {params.script}  -f {input} -o {output}

       """

rule bam_to_bed:
    input:
        "results/2.mapped/{sample}.sorted.bam"
    output:
        "results/5.bed/{sample}.bed"
    shell:
       """
        bedtools bamtobed -split -i {input} > {output}

       """

rule bed_to_cov:
    input:
        "results/5.bed/{sample}.bed"
        
    params:
        ref_bed="resources/refs/target_site.bed"
  
    output:
        "results/6.cov_per_site/{sample}.cov"
    shell:
       """
        bedtools coverage  -counts -a {params.ref_bed} -b {input} > {output}

       """
rule cov_tally:
    input:
        "results/6.cov_per_site/{sample}.cov"
    params:
        script = 'workflow/scripts/quick_cov_count.py'
    output:
        "results/7.tally/{sample}.tsv"
    shell:
       """
        $( which python ) {params.script}  -f {input} -o {output}

       """


























