#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio.Align import substitution_matrices
matrix_b90 = substitution_matrices.load("BLOSUM90")

import array
import re
import math
import sys

'''
Input file is a .csv file of the human aligned motif residue for each of the
    searched sites.
The code returns two PDF summary figures and a filtered csv file.
'''

match_code_dict = {"match":{"code":0, "label": "Match"},\
"swap": {"code":1, "label": "S/T Swap"}, "nopres": {"code":3, "label": "Not \n Phosphorylatable"},\
"gap": {"code": -1, "label": "Gap"}, "no_align": {"code":-2, "label": "No \n Residue"},\
"no_human": {"code": -3, "label": "No \n Reference"}}

def plot_cutoff(matches_only, mmatches_only, cutoff, cum_fdr, fig_FDR):

    width = 0.025
    fig, ax = plt.subplots()

    [m_counts,bins] = np.histogram(matches_only, bins= [i/20 for i in range(21)])
    [mm_counts,bins] = np.histogram(mmatches_only, bins= [i/20 for i in range(21)])

    ax.bar(x = bins[:len(bins)-1]-width/2.5, height = m_counts/sum(m_counts)*100,\
     width = 0.025, label='Matches (' + str(sum(m_counts)) + ")")
    ax.bar(x = bins[:len(bins)-1]+width/2.5, height = mm_counts/sum(mm_counts)*100,\
     width = 0.025, label='Mis-Matches (' + str(sum(mm_counts)) + ")")

    ax2 = ax.twinx()
    ax2.plot(bins[:len(bins)-1], np.array(cum_fdr)*100, color = "pink", label = "FDR", linewidth = '4')
    ax2.plot([cutoff + width/2 for i in range(25)],[i for i in range(25)], color = "black",\
         label = "Cutoff >= " + str(cutoff))

    ax2.set_ylabel('False Discovery Rate', weight = "bold")
    ax2.legend(bbox_to_anchor=(1.1, 0.4))
    ax2.set_ylim(bottom = 0)

    ax.set_ylabel("Percentage Scores (by class)", weight = "bold")
    ax.set_xlabel('Motif Score', weight = "bold")
    ax.set_title('Cutoff Determination Using Motif Scores \n of Matches and Mis-Matches', \
    weight = "bold")
    ax.legend(bbox_to_anchor=(1.1, 0.6))
    ax.set_ylim(bottom = 0)

    #plt.show()

    fig.savefig(fig_FDR, bbox_inches='tight')

def find_cutoff(matches_only, mmatches_only, FDR):

    bins_width = 0.1

    [m_counts,bins] = np.histogram(matches_only, bins= [i/20 for i in range(21)])
    [mm_counts,bins] = np.histogram(mmatches_only, bins= [i/20 for i in range(21)])

    cum_fdr = [sum(mm_counts[i:]) / (sum(mm_counts[i:])  + sum(m_counts[i:]))  for i in range(len(bins)-1)]

    cutoff = bins[np.argmax(np.array(cum_fdr) < FDR)]

    return cutoff, cum_fdr

def score_align(xen_in, hs_in, match_x, mismatch_x, matrix_b90):

    return_info = np.zeros((3,2))

    # ROWS
    # the first row is the alignment of Xenopus motif to itself
    # the second row is the alignment of the Human motif to itself
    # the last row is the alignment of Xenopus against Human

    # COLUMNS
    # the first column is the alignment going in the same direction
    # the second column is the alignemt going in the opposite direction

    ind_f = (0, 12)

    # assess both the forward and reverse alignments
    for i in range(len(ind_f)):

        # if going with aligned or reverse aligned
        which_way = ind_f[i]

        #for each position in the alignment except the phospho-residue
        for j in [0,1,2,3,4,5,7,8,9,10,11,12]:

            k = abs(j-which_way)

            # have different options to make sure that one of the aligned residues is not an "x" because
            # we are at the beginning or end of the protein

            # there are no "x"s
            if (xen_in[j] != "x" and xen_in[k] != "x" and hs_in[j] != "x" and hs_in[k] != "x"):

                return_info[2,i] += matrix_b90[xen_in[j]][hs_in[k]]
                return_info[0,i] += matrix_b90[xen_in[j]][xen_in[k]]
                return_info[1,i] += matrix_b90[hs_in[j]][hs_in[k]]

            # everything is an "x"
            elif (xen_in[j] == "x" and xen_in[k] == "x" and hs_in[j] == "x" and hs_in[k] == "x"):

                return_info[2,i] += match_x
                return_info[0,i] += match_x
                return_info[1,i] += match_x

            # if both Xenopus positions are "x," but neither human ones are
            elif ((xen_in[j] == "x" and xen_in[k] == "x") and (hs_in[j] != "x" and hs_in[k] != "x")):

                return_info[2,i] += mismatch_x
                return_info[0,i] += match_x
                return_info[1,i] += matrix_b90[hs_in[j]][hs_in[k]]

            # if both Human positions are "x," but neither xenopus one is
            elif ((xen_in[j] != "x" and xen_in[k] != "x") and (hs_in[j] == "x" and hs_in[k] == "x")):

                return_info[2,i] += mismatch_x
                return_info[0,i] += matrix_b90[xen_in[j]][xen_in[k]]
                return_info[1,i] += match_x

            # if one xenopus one is -- j != k
            elif (xen_in[j] == "x" or xen_in[k] == "x") or (hs_in[j] == "x" or hs_in[k] == "x"):

                if (xen_in[j] == "x" or xen_in[k] == "x") and (hs_in[j] == "x" or hs_in[k] == "x"):

                    return_info[2,i] += mismatch_x
                    return_info[0,i] += mismatch_x
                    return_info[1,i] += mismatch_x

                elif (xen_in[j] == "x" or xen_in[k] == "x"):

                    return_info[2,i] += mismatch_x
                    return_info[0,i] += mismatch_x
                    return_info[1,i] += matrix_b90[hs_in[j]][hs_in[k]]

                else:
                    return_info[2,i] += mismatch_x
                    return_info[0,i] += matrix_b90[xen_in[j]][xen_in[k]]
                    return_info[1,i] += mismatch_x

    return return_info

def plot_bar(match_summary, fig_bar, total_num):

    # make a bar plot of the match summary statistics
    fig, ax = plt.subplots()
    keys = list(match_summary.keys())
    ax.bar(keys, [match_summary[k]["freq"] for k in keys])
    ax.set(ylim = [0, 1])

    ax.set_ylabel("Fraction Total Residues")
    ax.set_title("Human Alignment Results for \n {} Phosphorylated Xenopus Residues".format(total_num))
    fig.savefig(fig_bar, bbox_inches='tight')


def score(match_file, fig_bar, FDR, fig_FDR, filtered_file):

    input_df = pd.read_csv(match_file, header = 0, low_memory = False, \
    dtype={"Human_Residue":str, "Match_Code":int})

    total_num = input_df.shape[0]

    match_summary  = {}
    for k,v in match_code_dict.items():
        match_summary[k] = {}
        match_summary[k]["count"] = sum(input_df["Match_Code"] == v["code"])
        match_summary[k]["freq"] = (match_summary[k]["count"])/total_num

    plot_bar(match_summary, fig_bar, total_num)

    # now find the match scores
    matched_motifs_F = input_df[input_df["Match_Code"]>=0]
    store_score = []

    # the value to use if both Xenopus and Human positions are "x"
    match_x = (np.diag(matrix_b90)).mean()
    store_mm_vals = []

    for i in range(len(matrix_b90)):
        for j in range(1,i):
            store_mm_vals.append(matrix_b90[i][j])

    # the value to use if one of the Xenopus or Human positions are "x"
    mismatch_x = np.mean(store_mm_vals)

    for i in range(len(matched_motifs_F)):

        xen_m = (matched_motifs_F.iloc[i,5])
        hs_m = (matched_motifs_F.iloc[i,6])

        if ("U" not in xen_m and "U" not in hs_m):

            return_info = score_align(xen_m, hs_m, \
            match_x, mismatch_x, matrix_b90)

            f_best = max(return_info[:,0])
            f_worst = min(return_info[:,1])
            store_score.append(round((return_info[2,0] - f_worst)/(f_best - f_worst),2))

        else:
            print("found Us")
            store_score.append(0)

    store_score_pd = pd.DataFrame(store_score)

    matched_motifs_F = matched_motifs_F.assign(Motif_Score = store_score_pd.values)

    # now we want to find the FDR cut off
    # the matches - perfect of S,T swap - match codes 0,1
    matches_only = matched_motifs_F[(matched_motifs_F["Match_Code"].isin([0, 2]))]

    # matches with alignment but no phosphorylateable residue
    mmatches_only = matched_motifs_F[matched_motifs_F["Match_Code"]==3]

    # find the cutoff
    cutoff, cum_fdr = find_cutoff(matches_only["Motif_Score"], mmatches_only["Motif_Score"], FDR)

    # make a plot
    plot_cutoff(matches_only["Motif_Score"], mmatches_only["Motif_Score"], cutoff, cum_fdr, fig_FDR)

    # return the matches that pass the threshold
    matched_motifs_pass = \
    matched_motifs_F.loc[(matched_motifs_F["Motif_Score"] >= cutoff) & \
    (matched_motifs_F["Match_Code"].isin([0,2])), :]

    matched_motifs_pass.to_csv(filtered_file, sep = ",", index = False)

if __name__ == '__main__':

    match_file = sys.argv[1]
    FDR = float(sys.argv[3])

    fig_bar = sys.argv[2]
    fig_FDR = sys.argv[4]

    filtered_file = sys.argv[5]

    score(match_file, fig_bar, FDR, fig_FDR, filtered_file)
