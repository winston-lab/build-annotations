#!/usr/bin/env python

import os

configfile: "config.yaml"

localrules:
    build_genic_annotation,
    build_convergent_annotation,
    build_divergent_annotation,
    build_intergenic_annotation,
    build_gc_coverage,
    build_motif_database,

#get all motif names from motif databases, cleaning nasty characters in some motif names
MOTIFS_DNA = set(subprocess.run(args="meme2meme " + " ".join(config["motifs"]["dna_motif_databases"]) + " | grep -e '^MOTIF' | cut -d ' ' -f2 | sed 's/\//_/g; s/&/_/g; s/{/[/g; s/}/]/g' ", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()) if config["motifs"]["build_motif_databases"] and config["motifs"]["dna_motif_databases"] else ""
MOTIFS_RNA = set(subprocess.run(args="meme2meme " + " ".join(config["motifs"]["rna_motif_databases"]) + " | grep -e '^MOTIF' | cut -d ' ' -f2 | sed 's/\//_/g; s/&/_/g; s/{/[/g; s/}/]/g' ", shell=True, stdout=subprocess.PIPE, encoding='utf-8').stdout.split()) if config["motifs"]["build_motif_databases"] and config["motifs"]["rna_motif_databases"] else ""

wildcard_constraints:
    nucleotide="dna|rna"

rule all:
    input:
        "config.yaml",
        expand("annotations/" + config["genome"]["name"] + "_{category}-regions.bed", category=["genic", "convergent", "divergent", "intergenic"]),
        "gc_pct/" + os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0] + "-GC_pct.bw",
        "motifs/" + config["genome"]["name"] + "_all_dna_motifs.bed" if config["motifs"]["build_motif_databases"] and config["motifs"]["dna_motif_databases"] else [],
        "motifs/" + config["genome"]["name"] + "_all_rna_motifs.bed" if config["motifs"]["build_motif_databases"] and config["motifs"]["rna_motif_databases"] else []

rule build_genic_annotation:
    input:
        transcripts = config["genome"]["transcript_annotation"],
        orfs = config["genome"]["orf_annotation"],
        fasta = config["genome"]["fasta"]
    output:
        "annotations/" + config["genome"]["name"] + "_genic-regions.bed"
    params:
        windowsize = config["genic_windowsize"]
    conda:
        "envs/build_annotations.yaml"
    log:
        "logs/build_genic_annotation.log"
    shell: """
        (python scripts/build_genic_annotation.py -t {input.transcripts} -o {input.orfs} -d {params.windowsize} -g <(faidx {input.fasta} -i chromsizes) -p {output}) &> {log}
        """

rule build_convergent_annotation:
    input:
        transcripts = config["genome"]["transcript_annotation"],
    output:
        "annotations/" + config["genome"]["name"] + "_convergent-regions.bed"
    params:
        max_dist = config["max_convergent_dist"]
    conda:
        "envs/build_annotations.yaml"
    log:
        "logs/build_convergent_annotation.log"
    shell: """
        (awk -v adist={params.max_dist} 'BEGIN{{FS=OFS="\t"}} $6=="+" {{ if(($3-$2)>adist) print $1, $2, $2+adist, $4, $5, "-" ; else print $0 }} $6=="-" {{if (($3-$2)>adist) print $1, $3-adist, $3, $4, $5, "+"; else print $0}}' {input.transcripts} > {output}) &> {log}
        """

rule build_divergent_annotation:
    input:
        transcripts = config["genome"]["transcript_annotation"],
        fasta = config["genome"]["fasta"]
    output:
        "annotations/" + config["genome"]["name"] + "_divergent-regions.bed"
    params:
        max_dist = config["max_divergent_dist"]
    conda:
        "envs/build_annotations.yaml"
    log:
        "logs/build_divergent_annotation.log"
    shell: """
        (bedtools flank -l {params.max_dist} -r 0 -s -i {input.transcripts} -g <(faidx {input.fasta} -i chromsizes) | awk 'BEGIN{{FS=OFS="\t"}} $6=="+"{{print $1, $2, $3, $4, $5, "-"}} $6=="-"{{print $1, $2, $3, $4, $5, "+"}}' > {output}) &> {log}
        """

rule build_intergenic_annotation:
    input:
        transcripts = config["genome"]["transcript_annotation"],
        fasta = config["genome"]["fasta"]
    output:
        "annotations/" + config["genome"]["name"] + "_intergenic-regions.bed"
    params:
        genic_up = config["genic_windowsize"]
    conda:
        "envs/build_annotations.yaml"
    log:
        "logs/build_intergenic_annotation.log"
    shell: """
        (bedtools slop -s -l {params.genic_up} -r 0 -i {input.transcripts} -g <(faidx {input.fasta} -i chromsizes | sort -k1,1) | sort -k1,1 -k2,2n | bedtools complement -i stdin -g <(sort -k1,1 <(faidx {input.fasta} -i chromsizes)) > {output}) &> {log}
        """

rule build_gc_coverage:
    input:
        fasta = config["genome"]["fasta"],
    output:
        "gc_pct/" + os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0] + "-GC_pct.bw"
    params:
        binsize = config["gc_pct_window"] #must be odd integer
    conda:
        "envs/build_annotations.yaml"
    log:
        "logs/build_gc_coverage.log"
    shell: """
        (python scripts/gc_content.py -f {input.fasta} -w {params.binsize} -o {output}) &> {log}
        """

motif_databases = {"dna": config["motifs"]["dna_motif_databases"],
                   "rna": config["motifs"]["rna_motif_databases"]}

rule build_motif_database:
    input:
        fasta = config["genome"]["fasta"],
        motif_db = lambda wc: motif_databases.get(wc.nucleotide) if config["motifs"]["build_motif_databases"] and motif_databases.get(wc.nucleotide) else []
    output:
        "motifs/" + config["genome"]["name"] + "_all_{nucleotide}_motifs.meme"
    conda:
        "envs/build_annotations.yaml"
    log:
        "logs/build_motif_database_{nucleotide}.log"
    shell: """
        (meme2meme -bg <(fasta-get-markov {input.fasta}) {input.motif_db} | sed -e 's/\//_/g; s/&/_/g; s/{{/[/g; s/}}/]/g' > {output}) &> {log}
        """

#run fimo in parallel for each motif
rule fimo:
    input:
        fasta = config["genome"]["fasta"],
        motif_db = "motifs/" + config["genome"]["name"] + "_all_{nucleotide}_motifs.meme"
    output:
        bed = temp("motifs/.{nucleotide}_{motif}.bed") # a BED6+2 format
    params:
        alpha = config["motifs"]["fimo_pval"]
    conda:
        "envs/build_annotations.yaml"
    log:
        "logs/fimo/fimo_{nucleotide}_{motif}.log"
    shell: """
        (fimo --motif {wildcards.motif} --bgfile <(fasta-get-markov {input.fasta}) --thresh {params.alpha} --text {input.motif_db} {input.fasta} | awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $3, $4-1, $5, $1, -log($8)/log(10), $6, $2, $10}}' > {output.bed}) &> {log}
        """

rule cat_fimo_motifs:
    input:
        bed = lambda wc: expand("motifs/.{nucleotide}_{motif}.bed", nucleotide=wc.nucleotide, motif=(MOTIFS_RNA if wc.nucleotide=="rna" else MOTIFS_DNA))
    output:
        bed = "motifs/" + config["genome"]["name"] + "_all_{nucleotide}_motifs.bed"
    conda:
        "envs/build_annotations.yaml"
    log:
        "logs/cat_fimo_motifs_{nucleotide}.log"
    threads:
        config["threads"]
    shell: """
        (cat {input.bed} | sort -k1,1 -k2,2n --parallel={threads} > {output.bed}) &> {log}
        """

