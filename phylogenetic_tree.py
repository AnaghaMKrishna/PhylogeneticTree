#!/usr/bin/env python3

import upgma
import os

MAX_AMPLICON_SIZE = 2000

#path for primer and assembly directory
primer_path = os.path.join("./data/16s_tree/primers/", "general_16S_515f_806r.fna")
assembly_dir = "./data/16s_tree/assemblies/"

amplicon_list = []
amplicon_dict = {}

#Find amplicons using isPCR for each assembly file and store them in dictionary
for assembly_file in os.listdir(assembly_dir):
    assembly_path = os.path.join(assembly_dir, assembly_file)
    amplicon_list.extend(upgma.ispcr(primer_path, assembly_path, MAX_AMPLICON_SIZE).split("\n"))
    #skip over the header lines and add sequences into a set to remove duplicate amplicons
    amplicon_dict[assembly_file[:assembly_file.index(".")]] = list({seq for seq in amplicon_list[1::2]})
    amplicon_list.clear() #empty the list before next iteration

#uncomment to print amplicons for each organism
# for k,v in amplicon_dict.items():
#     print(f"{k}\n{v}\n\n")      

#Align amplicons using NW global alignment 
norm_sim_matrix, org_list = upgma.generate_dist_matrix(amplicon_dict)
#Find clusters based on their similarity score to create phylogenetic tree
newick_str = upgma.finding_clusters(norm_sim_matrix, org_list)
print(newick_str)