#!/usr/bin/env python

import os

configfile: "config.yaml"

localrules:
    build_genic_annotation,
    build_convergent_annotation,
    build_divergent_annotation,
    build_intergenic_annotation,
    build_gc_coverage
    build_motif_database,

rule all:
    input:
        "config.yaml",
        os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "genic-regions.bed",
        os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "convergent-regions.bed",
        os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "divergent-regions.bed",
        os.path.dirname(os.path.abspath(config["genome"]["transcripts"])) + "/" + config["genome"]["prefix"] + "intergenic-regions.bed",
        os.path.splitext(os.path.abspath(config["genome"]["fasta"]))[0] + "-GC_pct.bw"

rule build_genic_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        fasta = config["genome"]["fasta"]
    output:
        "annotations/" + config["genome"]["prefix"] + "genic-regions.bed"
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
        "annotations/" + config["genome"]["prefix"] + "convergent-regions.bed"
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
        "annotations/" + config["genome"]["prefix"] + "divergent-regions.bed"
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
        "annotations/" + config["genome"]["prefix"] + "intergenic-regions.bed"
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
        "gc_pct/" + os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0] + "-GC_pct.bw"
    params:
        binsize = config["gc-pct-window"] #must be odd integer
    conda: "envs/build_annotations.yaml"
    log: "logs/build_gc_coverage.log"
    shell: """
        (python scripts/gc_content.py -f {input.fasta} -w {params.binsize} -o {output}) &> {log}
        """

rule build_motif_database:
    input:
        fasta = config["genome"]["fasta"],
        motif_db = config["motifs"]["databases"]
    output:
        "motifs/" + config["genome"]["prefix"] + "allmotifs.meme"
    log: "logs/build_motif_database.log"
    shell: """
        (meme2meme -bg <(fasta-get-markov {input.fasta}) {input.motif_db} | sed -e 's/\//_/g; s/&/_/g; s/{{/[/g; s/}}/]/g' > {output}) &> {log}
        """

#run fimo in parallel for each motif
rule fimo:
    input:
        fasta = config["genome"]["fasta"],
        motif_db = "motifs/" + config["genome"]["prefix"] + "allmotifs.meme"
    output:
        bed = temp("motifs/.{motif}.bed") # a BED6+2 format
    params:
        alpha = config["motifs"]["fimo-pval"]
    log: "logs/fimo/fimo_{motif}.log"
    shell: """
        (fimo --motif {wildcards.motif} --bgfile <(fasta-get-markov {input.fasta}) --thresh {params.alpha} --text {input.motif_db} {input.fasta} | awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $3, $4-1, $5, $1, -log($8)/log(10), $6, $2, $10}}' > {output.bed}) &> {log}
        """

rule cat_fimo_motifs:
    input:
        bed = expand("motifs/.{motif}.bed", motif=MOTIFS)
    output:
        bed = "motifs/" + config["genome"]["prefix"] + "allmotifs.bed"
    log: "logs/cat_fimo_motifs.log"
    threads: config["threads"]
    shell: """
        (cat {input.bed} | sort -k1,1 -k2,2n --parallel={threads} > {output.bed}) &> {log}
        """
