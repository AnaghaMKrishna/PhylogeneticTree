#!/usr/bin/env python3

import upgma
import numpy as np

MATCH = 1
MISMATCH = -1
GAP = -1

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

def generate_dist_matrix(amplicon_dict: dict[list[str]]) -> tuple[np.ndarray, list[str]]:
    """
    Initialize and fill distance matrix with scores from pair-wise comparision of amplicons

    Args:
        amplicon_dict: Dictionary of amplicons

    Returns:
        2D distance matrix and list of organism names
    """
    #initialize similarity matrix with arbitrary big value
    dist_matrix = np.full((len(amplicon_dict), len(amplicon_dict)), 100.0)

    #list with names of all organisms
    org_list = list(amplicon_dict.keys())
    
    #align every organism's amplicon with every other to get the similarity score and fill the matrix
    for idx1, org1 in enumerate(amplicon_dict.keys()):
        for idx2, org2 in enumerate(amplicon_dict.keys()):
            if org1 == org2:
                continue
            else:
                #when both amplicons are empty str
                if len(amplicon_dict[org1]) == 0 and len(amplicon_dict[org2]) == 0:
                    score = 0.000
                    dissim_score = 1.000
                    dist_matrix[idx1, idx2] = dissim_score
        
                #when one of the amplicons are empty str
                elif len(amplicon_dict[org1]) == 0 or len(amplicon_dict[org2]) == 0:
                    score = len(amplicon_dict[org1])*GAP if len(amplicon_dict[org1]) > 0 else len(amplicon_dict[org2])*GAP
                    #invert similarity score to get dissimilarity score and normalize
                    dissim_score = max(len(seq1), len(seq2))*MATCH - score
                    dist_matrix[idx1, idx2] = np.round(dissim_score / max(len(seq1), len(seq2)), 3)
                    
                #when both amplicons are not empty
                else:
                    # seq1 = next(iter(amplicon_dict[org1]))
                    # seq2 = next(iter(amplicon_dict[org2]))
                    
                    #choose the first amplicon from the list
                    seq1 = amplicon_dict[org1][0]
                    seq2 = amplicon_dict[org2][0]

                    #NW global alignment
                    _, score1 = upgma.needleman_wunsch(seq1, seq2, MATCH, MISMATCH, GAP)
                    _, score2 = upgma.needleman_wunsch(seq1, reverse_complement(seq2), MATCH, MISMATCH, GAP)
                    #select max score 
                    score = score1 if score1 > score2 else score2
                    
                    # dist_matrix[idx1, idx2] = score
                    # dist_matrix[idx1, idx2] = np.round(score / max(len(seq1), len(seq2)), 3)

                    #invert similarity score to get dissimilarity score and normalize
                    dissim_score = max(len(seq1), len(seq2))*MATCH - score
                    dist_matrix[idx1, idx2] = np.round(dissim_score / max(len(seq1), len(seq2)), 3)

    #uncomment to print distance matrix                
    # for i in dist_matrix:
    #     print(i)

    return (dist_matrix, org_list)

def finding_clusters(dist_matrix: np.ndarray, org_list: list[str]) -> str:
    """
    Find related organisms using UPGMA

    Args:
        sim_matrix: 2D distance matrix
        org_list: List of organism names

    Returns:
        Newick format of phylogenetic tree of given organisms
    """
    #dictionary to store distance for constructing newick string
    org_dist = {}

    while len(dist_matrix) > 1:
        #find minimum distance and it's index
        idx_pair = np.argwhere(dist_matrix == dist_matrix.min())[0]
        idx1 = idx_pair[0]
        idx2 = idx_pair[1]
        min_value = dist_matrix[idx1][idx2]
        
        org1_avg = np.round((min_value / 2), 3)
        org2_avg = np.round((min_value / 2), 3)

        modified_dist = []
        #for row and col containing max value, calculate the average for every other row and col
        for i in range(0, len(dist_matrix)):
            if i == idx1 or i == idx2:
                continue
            else:
                modified_dist.append(np.round((dist_matrix[idx1][i] + dist_matrix[idx2][i]) / 2, 3))

        # print(modified_dist)

        #convert to numpy array and insert diagonal element
        np_modified_dist = np.array(modified_dist)
        np_modified_dist = np.insert(np_modified_dist, idx1, 100.0)

        #delete the row and col at idx2
        row_to_del = idx1
        col_to_del = idx2
        dist_matrix = np.delete(dist_matrix, np.s_[col_to_del], 0)
        dist_matrix = np.delete(dist_matrix, np.s_[col_to_del], 1)
        
        #replace values in idx1-th row and column with modified_dist
        dist_matrix[row_to_del, : ] = np_modified_dist[:]
        dist_matrix[:, row_to_del] = np_modified_dist.T[:]

        #calculate the edge length for the node to satisfy ultrametricity
        if org_list[idx1] in org_dist:
            org1_avg = np.round((min_value / 2) - org_dist[org_list[idx1]], 3)
        if org_list[idx2] in org_dist:
            org2_avg = np.round((min_value / 2) - org_dist[org_list[idx2]], 3)

        # print(dist_matrix)

        #generate newick string
        org_list[idx1] = '(' + org_list[idx1] + ":" + str(org1_avg) + "," + org_list[idx2] + ":" + str(org2_avg) + ')'
        org_dist[org_list[idx1]] = np.round((min_value) / 2, 3)
        del org_list[idx2]

    return org_list[0] + ";"

