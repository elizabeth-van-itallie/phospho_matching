#!/usr/bin/env python3
import pandas as pd
import numpy as np
import re
import math

def get_mod_res(to_regexp):

    # get just the residue number from "MOD_RSD"
    m = []
    out_res_num = []

    for i in range(len(to_regexp)):
        m = re.search('(?<=[T|S|Y])\w+',to_regexp[i])
        if (m is not None):
            out_res_num.append(m.group(0))
        else:
            #out_res_num.append(math.nan)
            out_res_num.append(0)

    return out_res_num

def get_Human_ID(ref_to_split):

    out_split = []
    for i in range(len(ref_to_split)): out_split.append(ref_to_split[i].split('|')[3])

    return np.array(out_split)

def get_info(pp_info_file, match_filter_file):

    # load the phospho site plus information as a data frame
    pp_info_PD = pd.read_csv(pp_info_file, sep='\t', skiprows = (0,1), header = 0)

    # get just the residue number
    #   and store it as a new column in the data frame
    print(pp_info_PD.head())
    out_res_num = get_mod_res(pp_info_PD.MOD_RSD.to_numpy())
    pp_info_PD["RES_NUMBER"] = np.array(out_res_num).astype("int64") # np.array(out_res_num)


    # load the match information
    match_info_PD = pd.read_csv(match_filter_file, sep = ",", header = 0)
    match_info_PD["Human_ACCID"] = get_Human_ID(match_info_PD.Human_Reference.to_numpy())
    # now we want to merge these dataframe on ACC-iD and residue numbers


    merged_PD = match_info_PD.merge(pp_info_PD, how = "left", \
    left_on = ["Human_ACCID", "Human_Residue"], \
    right_on = ["ACC_ID", "RES_NUMBER"])

    match_info_PD["PP_LT"] = merged_PD["LT_LIT"]
    match_info_PD["PP_AB"] = merged_PD["CST_CAT#"]

    match_info_PD.drop(columns = ["Human_ACCID"], inplace = True)

    return match_info_PD

if __name__ == "__main__":

    pp_info_file = snakemake.input[0]
    match_filter_file = snakemake.input[1]
    match_filter_info_file = snakemake.output[0]

    match_filter_info_pd = get_info(pp_info_file, match_filter_file)
    match_filter_info_pd.to_csv(match_filter_info_file)
