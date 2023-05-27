#!/usr/bin/env python3

import step3_match as st3

"""
For each alignment file parse the BTOP string to determine what
    (phospho-acceptor, not phospho-acceptor, gap, nothing) the phosphorylated
    xenopus residues align to
"""

import pathlib
import csv
import re
from Bio import SeqIO

btop_dict = {'QueryName': 0, 'SubjectName': 1, 'MatchPercent': 2, \
'AlignLength': 3, 'QueryStart': 6, 'QueryEnd': 7, 'SubjectStart': 8, \
'SubjectEnd': 9, 'BTOP': 11}

match_dict = {'Xen_Res': 'Xen_Res', 'Human_Res': 'Human_Res','Match': 'Match',\
'Human_Info': 'Human_Info'}

info_dict = {"file_del": "\t", "Xen_Reference": 1, "Xen_Residues": 2,\
 "Xen_Motifs": 3, "file_sep": ";"}

output_header = ["Xenopus_Reference","Human_Reference","Match_Code", \
 "Xen_Residue","Human_Residue","Xen_Motif","Human_Motif"]

# create a dictionary of motifs from a fasta file and residues
def make_motifs_dict(fasta_dict, phospho_dict):

    motif_dict = {}

    for seq_id, phos_refs in phospho_dict.items():
        motif_dict[seq_id] = {}

        for p_str in phos_refs:
            p = int(p_str)-1
            s_pad_mot = ''
            e_pad_mot = ''
            ref_len = len(fasta_dict[seq_id])

            start = max(p-6,0)
            if (start == 0):
                s_pad_mot = 'x'*(6-p)

            end = min(p+7,ref_len)
            if (end-p < 7):
                e_pad_mot = 'x'*(6-(end-p-1))

            # find the sequence part of the motif from the input reference
            mot_seq = fasta_dict[seq_id][start:end]

            # combine the AAs and the padding for final motif
            motif_dict[seq_id][p_str] = s_pad_mot + mot_seq + e_pad_mot

    return motif_dict


# read a file as a list of lists and skip the lines that start with #
def read_blastout_skipHASH(filename, separator, quote):
    """
    Inputs:
        filename  - name of file
        separator - character that separates fields
        quote     - character used to optionally quote fields
    Output:
        Returns a list of lists where each item in the list
    corresponds to a row in the file.
    """
    table = []
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile, delimiter=separator, quotechar=quote)
        for row in csvreader:
            # skip the lines that start with #
            if (row[0].startswith("#") == False):
                table.append(row)
    return table

# parse a blastp BTOP string to process the alignment results for each amino acid
def parse_btop(list_in, btop_dict):
    """
    Inputs:
    list_in - list with the line from BLASTP
    btop_dict - the dictionary with the indexes to the information in the blastp
        output lines
    Output:
    match_list - a list the size of a the query alignment length that codes what
        each residue of query aligns to in subject
        -2 - no alignment
        -1 - gap
        0 - match
        2 - S/T swap
        3 - mismatch
    sub_num - list that is the same size as match_list that has the residue #
        for the human sequence that corresponds to the match information
        ( if a residue in the subject aligns to a gap in the query that residues
        number will not be included in sub_num)
    """
    # store the coded match for alignment length
    match_list = []
    quer_num = []
    sub_num = []

    btop_orig = list_in[btop_dict["BTOP"]]
    ind_s = 0
    quer_s = int(list_in[btop_dict["QueryStart"]])
    sub_s = int(list_in[btop_dict["SubjectStart"]])

    while (ind_s < len(btop_orig)):
        # return a match object to see if there is a non-numeric character
        not_match = re.search(r"\D",btop_orig[ind_s:])
        ind_gap = btop_orig[ind_s:].find("-")

        if (not_match):

            # the next values are numbers -- MATCHES
            if (not_match.start(0) > 0):
                ind_e = ind_s + not_match.start(0)
                match_int = int(btop_orig[ind_s:ind_e])
                match_list.extend([0]*match_int)

                sub_num.extend([i for i in range(sub_s,sub_s+match_int)])

                ind_s = ind_e
                quer_s = quer_s + match_int
                sub_s = sub_s + match_int

            # if there is not a match need to determine if there is a GAP
            # if the first position is a gap then - SUBJECT gap
            elif (ind_gap == 0):
                # the next position is a SUBJECT gap
                #match_list.extend([1])
                #aa_list.extend(["0"])
                ind_s += 2
                sub_s += 1

            elif (ind_gap == 1):
            #     # the next position is a QUERY gap
                match_list.extend([-1])

                #quer_num.extend([quer_s])
                sub_num.extend([sub_s])

                ind_s += 2
                quer_s += 1

            else:
                if (re.match(r"[S|T][S|T]",btop_orig[ind_s:ind_s+2])):
                    match_list.extend([2])

                # the next position is a mismatch
                else:
                    match_list.extend([3])

                #quer_num.extend([quer_s])
                sub_num.extend([sub_s])

                ind_s += 2
                quer_s += 1
                sub_s += 1
        else:
            # this is the case where there are all matches at the end

            match_int = int(btop_orig[ind_s:])

            match_list.extend([0]*match_int)
            ind_s = len(btop_orig)

            sub_num.extend([i for i in range(sub_s,sub_s+match_int)])

    return match_list, sub_num # quer_num

# use the output of parse_btop to figure out (1) if xenopus residue aligned to a
# human residue and if so (2) what the number of that residue is
def check_match(match_list, sub_num, xen_phos_res, match_info, btop_dict):
    """
    Inputs:
        match_list: information about the alignment of the sequences
            -2 - no alignment
            -1 - gap
            0 - match
            2 - S/T swap
            3 - mismatch
        sub_num: the list of the human protein residues that are aligned
        (if there is alignment of a subject residue to a gap in query that residue number/s will be skipped)
        xen_phos_res: the list of the measured phosphorylated xenopus residues
        match_info: the BLASTP output line that has alignement details
        ( the query start and end amino acids are relevant for this code)
        btop_dict: the DICTIONARY that gives the indexes into match_info to access information
    Output:
        phos_match_list: the list that is the same size as xen_phos_res that
            has the codes from match_list about the alignment
        human_phos_res: list that is the same size as xen_phos_res that
            has the number to human reference for the phosphorylated residue
    """
    phos_match_list = []
    human_phos_res = []

    for ind in xen_phos_res:
        # determine if that index aligned
        if int(ind) in list(range(int(match_info[btop_dict["QueryStart"]]),\
        int(match_info[btop_dict["QueryEnd"]]))):
            # determine what the residue aligned too
            phos_match_list.append(match_list[int(ind) - int(match_info[btop_dict["QueryStart"]])])
            # determine what the human residue number is that it aligns to
            human_phos_res.append(sub_num[int(ind) - int(match_info[btop_dict["QueryStart"]])])

        else:
            # if the xenopus residue is not in the aligned sequence
            phos_match_list.append(-2)
            human_phos_res.append(-2)

    return phos_match_list, human_phos_res

# print the results of the alignment parsing into a text file
def print_align_results(match_info, btop_dict, xen_phos_res, phos_match_list, human_phos_res, output_file):

    with open(output_file, 'w', newline='') as csvfile:
        phos_match_writer = csv.writer(csvfile, delimiter=' ')

        phos_match_writer.writerow([match_info[btop_dict["QueryName"]], \
        match_info[btop_dict["SubjectName"]], match_info[btop_dict["MatchPercent"]], \
        match_info[btop_dict["AlignLength"]]])

        for ind in range(len(xen_phos_res)):
            phos_match_writer.writerow([xen_phos_res[ind], human_phos_res[ind], \
            phos_match_list[ind]])


def phospho_file_to_dict(phospho_file):

    phospho_dict = {}
    with open(phospho_file, mode = "r") as txt_file:
        site_reader = csv.reader(txt_file, delimiter = "\t")
        for info in site_reader:
            phospho_dict[info[0]] = sorted(list(info[1].split(";")))

    return phospho_dict

# for each BLASTP alignment
def for_each_match(phospho_dict, file_end_i, input_folder):
    """
    Inputs:
        table_in: the list of lists with information about the matched references
        qr_in: the index the sub-lists of table_in with the query name which can
            be used to access to create the blastp access file
    Output:
    """


    no_blast_alignment = [] # the names of the sequences that do not have blast output
    alignment_results = {} # dictonary with the information about alignment to non-match or gap

    for xen_ref in phospho_dict.keys():
        # determine the name of the blastp output file
        i_file_name = pathlib.Path(input_folder, xen_ref + file_end_i)
        out_text = read_blastout_skipHASH(i_file_name, "\t",'"')
        # will call another function for each of the lists in out_text
        if (len(out_text)>0):
            [match_list, sub_num] = parse_btop(out_text[0], btop_dict)


            # now find the entry in the xen_res dictionary
            #if row[qr_ind] in xen_res.keys():
            xen_phos_res = phospho_dict[xen_ref]

            [phos_match_list, human_phos_res] = \
            check_match(match_list, sub_num, xen_phos_res, \
            out_text[0], btop_dict)

            sub_dict = {}
            sub_dict[match_dict["Xen_Res"]] = xen_phos_res
            sub_dict[match_dict["Human_Res"]] = human_phos_res
            sub_dict[match_dict["Match"]] = phos_match_list
            sub_dict[match_dict["Human_Info"]] = \
            out_text[0][btop_dict["SubjectName"]]

            alignment_results[xen_ref] = sub_dict

        else:
            no_blast_alignment.append(xen_ref)

    return no_blast_alignment, alignment_results

if __name__ == "__main__":

    phospho_sites_file = snakemake.input[0]
    input_folder = snakemake.input[1]
    file_end_i = "-out.txt"

    phospho_dict = phospho_file_to_dict(phospho_sites_file)

    [no_blast_alignment, alignment_results] = \
    for_each_match(phospho_dict, file_end_i, input_folder)

    fasta_in = snakemake.input[2]
    fasta_in_dict = {rec.id: str(rec.seq)\
     for rec in SeqIO.parse(fasta_in, "fasta")}

    # get the motifs
    motifs_dict = make_motifs_dict(fasta_in_dict, phospho_dict)

    fasta_human_no_iso = snakemake.input[3]
    no_iso_dict = {rec.id: str(rec.seq)\
     for rec in SeqIO.parse(fasta_human_no_iso, "fasta")}

    output_file = snakemake.output[0]

    alignment_results_more = st3.compile_results(alignment_results, no_blast_alignment,\
    match_dict, motifs_dict, info_dict, no_iso_dict)

    st3.writeoutput_table(alignment_results_more, output_file, output_header)
