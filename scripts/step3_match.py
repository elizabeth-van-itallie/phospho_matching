#!/usr/bin/env python3

import csv
import re

def compile_results(alignment_results, no_alignment, match_dict, motif_dict,\
 xen_info_dict, human_fasta_dict):
    """
    Create a dictionary with keys that are xenopus references that has all the
    information about the alignment results with motifs where relevant.
    """

    rep_stop = re.compile("\*")

    # create a dictionary of the xenopus residues and motifs for each reference
    #xen_info_list = read_xeninfo_as_dict(xen_info_file, xen_info_dict["file_del"],\
    # '"', xen_info_dict["Xen_Reference"], xen_info_dict["Xen_Residues"],\
    # xen_info_dict["Xen_Motifs"], xen_info_dict["file_sep"])

    # for all of the xenopus residues that align to a human residue find the motif for the
    # human residue

    for xen_ref in alignment_results.keys():

        match_result = alignment_results[xen_ref][match_dict["Match"]]

        human_ref = alignment_results[xen_ref][match_dict["Human_Info"]]
        human_motif = []
        xen_motif = []
        xen_residues = alignment_results[xen_ref][match_dict["Xen_Res"]]

        #s_pad_mot = ''
        #e_pad_mot = ''

        for match_ind in range(len(match_result)):

            if (match_result[match_ind]>=0):
            # if there was alignment to a residue find the +/- 6 motif
                s_pad_mot = ''
                e_pad_mot = ''
            # FIRST construct the indexing into the human motif
                ref_len = len(human_fasta_dict[human_ref])
                align_res = alignment_results[xen_ref][match_dict["Human_Res"]][match_ind]
                # substrat one because the indexing for alignment started at 1
                align_ind = int(align_res)-1

                # need to check if the residue is so near to the N or C terminus
                #   that the  motif needs to be padded with "x"s to have full length
                start = max(align_ind-6,0)
                if (start == 0):
                    s_pad_mot = 'x'*(6-align_ind)

                end = min(align_ind+7,ref_len)
                if (end-align_ind < 7):
                    e_pad_mot = 'x'*(6-(end-align_ind-1))

                # find the sequence part of the motif from the human reference
                mot_seq = human_fasta_dict[human_ref][start:end]

                # combine the AAs and the padding for final motif
                mot_all = s_pad_mot + mot_seq + e_pad_mot

                # save the human motif
                human_motif.append(mot_all)

                # get the human motif from the xenopus dictionary and save it
                xen_motif.append(\
                motif_dict[xen_ref][xen_residues[match_ind]].replace("*", "x"))
                #xen_motif.append(xen_info_list[xen_ref][xen_residues[match_ind]])

            # if the xenopus residue doesn't align to a human res than add an
            #   empty entry in to the list of motifs
            else:
                human_motif.append('')
                xen_motif.append('')

        # after all the xenopus residues have been stepped through save the
        #   xenopus and human motifs in the dictionary
        alignment_results[xen_ref]["Human_Motifs"] = human_motif
        alignment_results[xen_ref]["Xen_Motifs"] = xen_motif


    for xen_ref in no_alignment:
        x_motifs = []
        x_res = []
        score_no = []

        in_dict = motif_dict[xen_ref]
        alignment_results[xen_ref] = {}
        alignment_results[xen_ref]["Xen_Res"] = list(in_dict.keys())
        alignment_results[xen_ref]["Match"] = len(in_dict.keys())*[-3]

    return alignment_results

def writeoutput_table(alignment_results_in, output_file, output_header):

    with open(output_file, 'w', newline='') as csvfile:
            out_writer = csv.writer(csvfile, delimiter=',')
            out_writer.writerow(output_header)

            for xen_ref in alignment_results_in.keys():

                if max(alignment_results_in[xen_ref]["Match"])>-2:

                    for i in range(len(alignment_results_in[xen_ref]["Xen_Res"])):

                        out_writer.writerow([xen_ref,\
                        alignment_results_in[xen_ref]["Human_Info"],\
                        alignment_results_in[xen_ref]["Match"][i],\
                    alignment_results_in[xen_ref]["Xen_Res"][i],\
                    str(alignment_results_in[xen_ref]["Human_Res"][i]),\
                    alignment_results_in[xen_ref]["Xen_Motifs"][i],\
                    alignment_results_in[xen_ref]["Human_Motifs"][i]])

                else:
                       for i in range(len(alignment_results_in[xen_ref]["Xen_Res"])):

                            out_writer.writerow([xen_ref,"",\
                            alignment_results_in[xen_ref]["Match"][i],\
                            alignment_results_in[xen_ref]["Xen_Res"][i],"",\
                            "",""])
