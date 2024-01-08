#!/usr/bin/env python3

import upgma
import os
import numpy as np

MAX_AMPLICON_SIZE = 2000

#directory path for primer and assembly files
primer_dir = "./data/Pseudomonas_tree/primers/"
assembly_dir = "./data/Pseudomonas_tree/assemblies/"

amplicon_list = []
primer_dict = {}
org_dict = {}

#iterate over each organism and primer file to find the amplicons and store them in nested dictionary
for primer_file in os.listdir(primer_dir):
    primer_path = os.path.join(primer_dir, primer_file)
    primer_name = primer_file[:primer_file.index(".")]
    primer_dict[primer_name] = {}

    for assembly_file in os.listdir(assembly_dir):
        assembly_path = os.path.join(assembly_dir, assembly_file)
        org_name = assembly_file[:assembly_file.index(".")]
        amplicon_list.extend(upgma.ispcr(primer_path, assembly_path, MAX_AMPLICON_SIZE).split("\n"))
        #skip over the header lines and add sequences into a set to remove duplicate amplicons
        org_dict[org_name] = list({amp for amp in amplicon_list[1::2]})
        amplicon_list.clear() #clear list before next iteration

    #for each primer, add primer name as key and value as dict containing organism and resp. amplicons
    primer_dict[primer_name] = org_dict.copy()
    org_dict.clear() #clear dict before next iteration

#uncomment to print amplicons for every organism for each primer
# for k,v in primer_dict.items():
#     print(f"{k}\n{v}\n\n")

#2D NumPy array to store the sum of similarity scores for each organism's amplicon for each primer
conc_dist_mat = np.zeros((len(primer_dict[primer_name]), len(primer_dict[primer_name])))

for org_ampl in primer_dict.values():
    dist_mat, org_list = upgma.generate_dist_matrix(org_ampl)
    conc_dist_mat = np.add(conc_dist_mat, dist_mat)

#uncomment to print the distance matix for concatenated amplicons
# print(conc_dist_mat)

#perform UPGMA to generate phylogenetic tree in newick format
phy_tree = upgma.finding_clusters(conc_dist_mat, org_list)

print(phy_tree)