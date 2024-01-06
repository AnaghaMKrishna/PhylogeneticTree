#!/usr/bin/env python3

import upgma
import os
# from itertools import combinations
import numpy as np

max_amplicon_size = 2000
match = 1
mismatch = -1
gap = -1

def reverse_complement(string:str) -> str:
    """
    Generate the reverse complement of the string using complement_dict

    Args:
        string : string to be reverse complemented

    Returns:
        reverse complemented string
    """
    complement_dict = {
    "A":"T",
    "G":"C",
    "T":"A",
    "C":"G"
    }
    return "".join(complement_dict.get(base) for base in string[::-1])

def call_ispcr() -> dict[set]:
    #needs to take path from user as args
    amplicon_list = []
    amplicon_dict = {}

    #retrieve the path of the primer file
    primer_path = os.path.join("./data/16s_tree/primers/", "general_16S_515f_806r.fna")

    #retrieve the path of the assembly files to perform isPCR on each assembly
    assembly_dir = "./data/16s_tree/assemblies/"
    for assembly_file in os.listdir(assembly_dir):
        assembly_path = os.path.join(assembly_dir, assembly_file)
        amplicon_list.extend(upgma.ispcr(primer_path, assembly_path, max_amplicon_size).split("\n"))
        #skip over the header lines and add sequences into a set to remove duplicate amplicons
        amplicon_dict[assembly_file[:assembly_file.index(".")]] = list({seq for seq in amplicon_list[1::2]})
        #set does not preserve the order of the amplicons, so each time NW is run, it gives different score
        #trying with list
        amplicon_dict[assembly_file[:assembly_file.index(".")]] = [seq for seq in amplicon_list[1::2]]
        amplicon_list.clear() #empty the list before next iteration

    #uncomment to print amplicons for each organism
    # for k,v in amplicon_dict.items():
    #     print(f"{k}\n{v}\n\n")      

    return amplicon_dict  

def generate_sim_matrix(amplicon_dict: dict[set]) -> list[list[int]]:
    #TD:add decorator to deal with multiple amplicons for same organism
    #Initialize similarity matrix with 0s
    sim_matrix = [[0 for _ in range(len(amplicon_dict))] for _ in range(len(amplicon_dict))]

    # print(type(amplicon_dict.keys()))
    #align every organism's amplicon with every other to get the similarity score and fill the matrix
    for ind1, org1 in enumerate(amplicon_dict.keys()):
        for ind2, org2 in enumerate(amplicon_dict.keys()):
            if org1 == org2:
                continue
            else:
                # seq1 = next(iter(amplicon_dict[org1]))
                # seq2 = next(iter(amplicon_dict[org2]))
                seq1 = amplicon_dict[org1][0]
                seq2 = amplicon_dict[org2][0]
                # print(f"{org1}: {seq1}")
                # print(f"{org2}: {seq2}")
                aln1, score1 = upgma.needleman_wunsch(seq1, seq2, match, mismatch, gap)
                aln2, score2 = upgma.needleman_wunsch(seq1, reverse_complement(seq2), match, mismatch, gap)
                # print(score)
                score = score1 if score1 > score2 else score2
                # sim_matrix[ind1][ind2] = score
                #normalize the scores
                sim_matrix[ind1][ind2] = round(score / max(len(seq1), len(seq2)), 3)
                # sim_matrix[ind1][ind2] = max(len(seq1), len(seq2))*match - score
                # break
        # break
    # for i in sim_matrix:
    #     print(i)
    return sim_matrix

def normalize_score(sim_matrix: list[list[float]]) ->list[list[float]]:
    #TD: calculate proper normal score
    for ind1, val in enumerate(sim_matrix):
        for ind2, ele in enumerate(val):
            sim_matrix[ind1][ind2] = sim_matrix[ind1][ind2]/100
    for i in sim_matrix:
        print(i)
    print(repr(sim_matrix))
    return sim_matrix

def finding_clusters(sim_matrix: list[list[float]]):
    # sim_matrix = [[0.0, 1.48, 1.58, 0.0, 0.03, 0.02, 1.59, 0.72, 0.08, 0.09], [1.48, 0.0, 1.53, 0.11, 0.13, 0.17, 1.3, 0.72, 0.14, 0.14], [1.58, 1.53, 0.0, -0.06, 0.08, 0.08, 1.58, 0.83, 0.11, 0.04], [0.0, 0.11, -0.06, 0.0, 0.77, 0.91, -0.19, -0.06, 0.68, 1.56], [0.03, 0.13, 0.08, 0.77, 0.0, 2.16, 0.04, 0.19, 1.28, 0.75], [0.02, 0.17, 0.08, 0.91, 2.16, 0.0, 0.04, 0.16, 1.48, 0.84], [1.59, 1.3, 1.58, -0.19, 0.04, 0.04, 0.0, 0.89, 0.12, 0.01], [0.72, 0.72, 0.83, -0.06, 0.19, 0.16, 0.89, 0.0, 0.14, -0.04], [0.08, 0.14, 0.11, 0.68, 1.28, 1.48, 0.12, 0.14, 0.0, 0.8], [0.09, 0.14, 0.04, 1.56, 0.75, 0.84, 0.01, -0.04, 0.8, 0.0]]
    org_list = ["Pseudomonas_aeruginosa_UCBPP-PA14", "Mycobacterium_tuberculosis_H37Rv", "Staphylococcus_aureus_NCTC_8325", "Sulfolobus_islandicus_M", "Vibrio_cholerae_N16961", "Escherichia_coli_K12", "Wolbachia", "Ferroplasma_acidiphilum_Y", "Treponema_phagedenis_strain_B43", "Nitrososphaera_viennensis_EN76"]
    np_sim_matrix = np.array(sim_matrix)
    for i in sim_matrix:
        print(i)
    # print(max_ind[0])
    # print(np_sim_matrix[max_ind[0,0],max_ind[0,1]])
    while len(np_sim_matrix) > 1:
        print(np_sim_matrix)
        max_indices = np.argwhere(np_sim_matrix == np_sim_matrix.max())
        max_value = np_sim_matrix[max_indices[0,0]][max_indices[0,1]]
        modified_sim = [] #np.empty(len(np_sim_matrix)-1)
        print(max_indices)
        print(max_value)
        # for the row and col containing max value, calculate the average for every other row and col
        for max_ind_pair in max_indices[::2]:
            max_avg = max_value / 2
            print(max_ind_pair)
            for i in range(len(np_sim_matrix)):
                if i == max_ind_pair[0] or i == max_ind_pair[1]:
                    continue
                else:
                    modified_sim.append(round((np_sim_matrix[max_ind_pair[0]][i-1] + np_sim_matrix[max_ind_pair[1]][i-1]) / 2, 3))
            # modified_sim[max_ind_pair[0]] = 0.00
            # print(modified_sim)
        # for row in range(len(np_sim_matrix)-1):
        #     for col in range(len(np_sim_matrix)):
        # reduced_sim_matrix = np.zeros((len(np_sim_matrix)-1, len(np_sim_matrix)-1))
            np_modified_sim = np.array(modified_sim)
            np_modified_sim = np.insert(np_modified_sim, max_ind_pair[0], 0.0)
            print(np_modified_sim)
            
            # np_mod_sim = np.insert(np_modified_sim, max_ind_pair[0], 0.0)
            # print(np_mod_sim)
            row_col_to_del = max_ind_pair
            # np_sim_matrix = np.delete(np_sim_matrix, np.s_[row_col_to_del[0]:row_col_to_del[1]:], 0)
            # np_sim_matrix = np.delete(np_sim_matrix, np.s_[row_col_to_del[0]:row_col_to_del[1]:], 1)
            np_sim_matrix = np.delete(np_sim_matrix, np.s_[row_col_to_del[1]], 0)
            np_sim_matrix = np.delete(np_sim_matrix, np.s_[row_col_to_del[1]], 1)
            # np_sim_matrix1 = np.insert(np_sim_matrix, np.s_[row_col_to_del[0]], np_mod_sim, axis = 0)
            # np_sim_matrix2 = np.insert(np_sim_matrix1, np.s_[row_col_to_del[0]], np_mod_sim, axis = 1)
            # # np_sim_matrix[:, row_col_to_del[0]:row_col_to_del[1] + 1]
            print(np_sim_matrix)
            
            #replace values in max index[0] with np_mod_sim
            np_sim_matrix[row_col_to_del[0], : ] = np_modified_sim[:]
            np_sim_matrix[:, row_col_to_del[0]] = np_modified_sim.T[:]
            # np_sim_matrix = np.insert(np_sim_matrix, row_col_to_del[0], np_modified_sim[:], axis = 0)
            # np_sim_matrix = np.insert(np_sim_matrix.T, row_col_to_del[0], np_modified_sim[:], axis = 0)


            # np_sim_matrix[:, row_col_to_del[0]:row_col_to_del[1] + 1] = np_modified_sim[:-1, np.newaxis]
            print(np_sim_matrix)

            #generate newick string
            org_list[max_ind_pair[0]] = '(' + org_list[max_ind_pair[0]] + ":" + str(max_avg) + "," + org_list[max_ind_pair[1]] + ":" + str(max_avg) + ')'
            del org_list[max_ind_pair[1]]
            print(org_list)
        



if __name__ == "__main__":
    amplicon_dict = call_ispcr()
    norm_sim_matrix = generate_sim_matrix(amplicon_dict)
    # norm_sim_matrix = normalize_score(sim_matrix)
    finding_clusters(norm_sim_matrix)




#check the optimal alignment for amplicons and its reverse complements
# aln1, score1 = upgma.needleman_wunsch(amplicon1, amplicon2, match, mismatch, gap)
# aln2, score2 = upgma.needleman_wunsch(reverse_complement(amplicon1), amplicon2, match, mismatch, gap)
# aln3, score3 = upgma.needleman_wunsch(amplicon1, reverse_complement(amplicon2), match, mismatch, gap)
# aln4, score4 = upgma.needleman_wunsch(reverse_complement(amplicon1), reverse_complement(amplicon2), match, mismatch, gap)

#check for maximum score and print the aligned sequences and corresponding score
# max_score = [[score1, aln1], [score2, aln2], [score3, aln3], [score4, aln4]]
# opt_aln = max_score.index(max(max_score))
# print("\n".join(max_score[opt_aln][1]))
# print(max_score[opt_aln][0])
