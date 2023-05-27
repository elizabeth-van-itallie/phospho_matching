#!/usr/bin/env python3

"""
    Code to generate the files for matching,
    call BLASTP for the paired Xenopus Human matches via the command line,
    and parse the file with Xenopus phospho-residues as dictionary
"""

import csv
import fastaparser
from pathlib import Path
import pathlib


# generate the single reference fasta files for the matched xenopus references
def generate_files4blast_xen():
    """
    Inputs:
        dict_in - the dictionary where the keys are the xenopus reference names that have phosphrylated residue
        fasta_in - the fasta file that has the xenopus reference information
        path_for_files
    Outputs:
        Nothing is returned
    """

    fasta_in = snakemake.input[0]
    phospho_sites = snakemake.input[1]
    path_for_files = snakemake.input[2]

    # read the fasta file into a dictionary
    phospho_dict = {}
    with open(phospho_sites, mode = "r") as txt_file:
        site_reader = csv.reader(txt_file, delimiter = "\t")
        for info in site_reader:
            phospho_dict[info[0]] = sorted(list(info[1].split(";")))

    fasta_dict = {}
    with open(fasta_in) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
        # seq is a FastaSequence object
            fasta_dict[seq.id] = seq.sequence_as_string()

    for key in phospho_dict.keys():
        if key in fasta_dict.keys():
        # generate the file
        #file_name_1 = row[col_ref] + ".fa"
            file_name = pathlib.Path(path_for_files, key + ".fa")
        # use fastaparser to write the single entry
            with open(file_name, 'w') as fasta_file:
                writer = fastaparser.Writer(fasta_file)
                writer.writefasta((key,fasta_dict[key]))

    with open(snakemake.output[0], "w") as f:
        f.write("single fastas are written")

if __name__ == '__main__':
    generate_files4blast_xen()
