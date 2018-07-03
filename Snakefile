#!/usr/bin/env python

import os

configfile: "config.yaml"

localrules:
    build_genic_annotation,
    build_convergent_annotation,
    build_divergent_annotation,
    build_intergenic_annotation,
    build_gc_coverage

rule all:
    input:
        genic_regions = os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "genic-regions.bed",
        convergent_regions = os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "convergent-regions.bed",
        divergent_regions = os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "divergent-regions.bed",
        intergenic_regions = os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "intergenic-regions.bed"

rule build_genic_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        fasta = config["genome"]["fasta"]
    output:
        os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "genic-regions.bed"
    params:
        windowsize = config["genic-windowsize"]
    conda: "envs/build_annotations.yaml"
    log : "logs/build_genic_annotation.log"
    shell: """
        (python scripts/build_genic_annotation.py -t {input.transcripts} -o {input.orfs} -d {params.windowsize} -g <(faidx {input.fasta} -i chromsizes) -p {output}) &> {log}
        """

rule build_convergent_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
    output:
        os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "convergent-regions.bed"
    params:
        max_dist = config["max-convergent-dist"]
    conda: "envs/build_annotations.yaml"
    log: "logs/build_convergent_annotation.log"
    shell: """
        (awk -v adist={params.max_dist} 'BEGIN{{FS=OFS="\t"}} $6=="+" {{ if(($3-$2)>adist) print $1, $2, $2+adist, $4, $5, "-" ; else print $0 }} $6=="-" {{if (($3-$2)>adist) print $1, $3-adist, $3, $4, $5, "+"; else print $0}}' {input.transcripts} > {output}) &> {log}
        """

rule build_divergent_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        fasta = config["genome"]["fasta"]
    output:
        os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "divergent-regions.bed"
    params:
        max_dist = config["max-divergent-dist"]
    conda: "envs/build_annotations.yaml"
    log: "logs/build_divergent_annotation.log"
    shell: """
        (bedtools flank -l {params.max_dist} -r 0 -s -i {input.transcripts} -g <(faidx {input.fasta} -i chromsizes) | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $3, $4, $5, "-"}} $6=="-"{{print $1, $2, $3, $4, $5, "+"}}' > {output}) &> {log}
        """

rule build_intergenic_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        fasta = config["genome"]["fasta"]
    output:
        os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "intergenic-regions.bed"
    params:
        genic_up = config["genic-windowsize"]
    conda: "envs/build_annotations.yaml"
    log: "logs/build_intergenic_annotation.log"
    shell: """
        (bedtools slop -s -l {params.genic_up} -r 0 -i {input.transcripts} -g <(faidx {input.fasta} -i chromsizes | sort -k1,1) | sort -k1,1 -k2,2n | bedtools complement -i stdin -g <(sort -k1,1 <(faidx {input.fasta} -i chromsizes)) > {output}) &> {log}
        """

rule build_gc_coverage:
    input:
        fasta = config["genome"]["fasta"],
    output:
        os.path.splitext(os.path.abspath(config["genome"]["fasta"]))[0] + "-GC_pct.bw"
    params:
        binsize = 11 #must be odd integer
    conda: "envs/build_annotations.yaml"
    log: "logs/build_gc_coverage.log"
    shell: """
        python scripts/gc_content.py -f {input.fasta} -w {params.binsize} -o {output}
        """

