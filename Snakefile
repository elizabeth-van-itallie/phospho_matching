#!/usr/bin/env python3

# export PATH=$PATH:/Users/evanitallie/ncbi-blast-2.5.0+/bin/
# conda activate snakemake
# snakemake --cores 3

import configparser
import subprocess

config = configparser.ConfigParser()
config.read("phos_match.config")

file_sites_in = config.get("phos_match_config_REAL", "file_sites_in")
fasta_file = config.get("phos_match_config_REAL", "fasta_file")
human_fasta_noIso = config.get("phos_match_config_REAL", "human_fasta_noIso")

fdr_cutoff = float(config.get("phos_match_config_REAL", "fdr_cutoff"))
b_eval = float(config.get("phos_match_config_REAL", "b_eval"))
number_workers = int(config.get("phos_match_config_REAL", "number_workers"))

output_match_file = config.get("phos_match_config_REAL", "output_match_file")
filtered_match_file = config.get("phos_match_config_REAL", "filtered_match_file")
pp_info_file = config.get("phos_match_config_REAL", "pp_info_file")
filtered_match_info_file = config.get("phos_match_config_REAL",\
 "filtered_match_info_file")

blast_output_files = "blast_output_files"
blast_input = "blast_input_files"

rule all:
    input:
        filtered_match_info_file,
        "checkpoint_files/packages_loaded.txt"
    shell:
        "echo snakemake pipeline done"

rule setup_folders:
    input:
        "checkpoint_files/packages_loaded.txt"
    output:
        directory(blast_input),
        directory(blast_output_files)
    params:
        folders = [blast_output_files, blast_input, "outputs", "checkpoint_files"]
    script:
        "scripts/setup_folders.py"

rule install_packages:
    output:
        "checkpoint_files/packages_loaded.txt"
    script:
        "scripts/install_packages.py"

rule human_no_iso:
    input:
        "fastas/human-phosphosite-fastas.fasta"
    output:
        human_fasta_noIso
    script:
        "scripts/to_no_iso.py"

rule make_single_fastas:
    input:
        fasta_file,
        file_sites_in,
        blast_input
    output:
        "checkpoint_files/single_fastas.txt"
    script:
        "scripts/blast_file_setup.py"

rule blast:
    input:
        "checkpoint_files/single_fastas.txt"
    output:
        "checkpoint_files/blast.txt"
    params:
        input_phos = file_sites_in,
        input_phos_ref_col = 0,
        input_phos_res_col = 1,
        input_phos_sep = "'\t'",
        input_files = blast_input,
        human_fasta = human_fasta_noIso,
        output_files = blast_output_files,
        eval = b_eval,
        blast_FMT = "'7 qacc sacc pident length mismatch gapopen qstart qend sstart send evalue btop'",
        num_workers = number_workers

    shell:
        "python3 scripts/step1_blast.py {params.input_phos} \
         {params.input_phos_ref_col} {params.input_phos_res_col} \
          {params.input_phos_sep} {params.input_files} {params.human_fasta}\
           {params.output_files} {params.eval} \
            {params.blast_FMT} {params.num_workers}"

rule compile_results:
     input:
        file_sites_in,
        blast_output_files,
        fasta_file,
        human_fasta_noIso,
        "checkpoint_files/blast.txt"
     output:
        output_match_file
     script:
        "scripts/step2_match.py"

rule motif_score:
    input:
        output_match_file
    params:
        input_file = output_match_file,
        FDR = fdr_cutoff,
        output_bar = "outputs/sites_bar.pdf",
        output_FDR = "outputs/cutoff_FDR.pdf",
        output_file = filtered_match_file
    output:
        "outputs/sites_bar.pdf",
        "outputs/cutoff_FDR.pdf",
        filtered_match_file
    shell:
        "python3 scripts/motif_score.py {params.input_file} {params.output_bar} \
         {params.FDR} {params.output_FDR} {params.output_file}"

rule get_pp_info:
    input:
        pp_info_file,
        filtered_match_file
    output:
        filtered_match_info_file
    script:
        "scripts/get_pp_info.py"

# have a clean RULE 
