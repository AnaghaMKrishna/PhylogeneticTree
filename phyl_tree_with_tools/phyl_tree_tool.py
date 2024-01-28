#!/usr/bin/env python3

import subprocess
import os

PCT_MATCH = 70.00
MAX_AMPLICON_SIZE = 2000

#directory path for primer and assembly files
#for 16S primer
primer_path = os.path.join("./data/16s_tree/primers/", "general_16S_515f_806r.fna")
assembly_dir = "./data/16s_tree/assemblies/"

#for multiple primers
primer_dir_multi = "./data/Pseudomonas_tree/primers/"
assembly_dir_multi = "./data/Pseudomonas_tree/assemblies/"

def single_locus(primer_path: str, assembly_dir:str) -> str:
    """
    generate phylogenetic tree for given assemblies considering single locus

    Args:
        primer_path: path to primer file
        assembly_dir: path to directory containing assemblies

    Returns:
        tree in newick format
    """
    amplicon_list = []
    amplicon_dict = {}

    #find amplicons using isPCR for each assembly file and store them in dictionary
    for assembly_file in os.listdir(assembly_dir):
        assembly_path = os.path.join(assembly_dir, assembly_file)
        org_name = assembly_file[:assembly_file.index(".")]
        amplicon_list.extend(ispcr(primer_path, assembly_path, MAX_AMPLICON_SIZE).split("\n"))
        #skip over the header lines and add sequences into a set to remove duplicate amplicons
        amplicon_dict[org_name] = list({seq for seq in amplicon_list[1::2]})
        amplicon_list.clear() #empty the list before next iteration

    #write amplicons to a file
    with open("amplicons.fna", "a") as fout:
            for org,amp in amplicon_dict.items():
                seq = f">{org}\n{amp[0]}\n"
                fout.write(seq)

    #perform MSA and generate newick tree
    aligned_file = msa("amplicons.fna", "aligned_ampl.efa")
    newick_str = newick_tree(aligned_file, "newick_tree.txt")

    return newick_str

def multi_loci(primer_dir_multi: str, assembly_dir_multi:str) -> str:
    """
    generate phylogenetic tree for given assemblies considering multiple loci

    Args:
        primer_path: path to directory containing primers
        assembly_dir: path to directory containing assemblies

    Returns:
        tree in newick format
    """
    amplicon_list = []
    primer_dict = {}
    org_dict = {}
    primer_list = []

    #iterate over each organism and primer file to find the amplicons and store them in nested dictionary
    for primer_file in os.listdir(primer_dir_multi):
        primer_path = os.path.join(primer_dir_multi, primer_file)
        primer_name = primer_file[:primer_file.index(".")]
        primer_list.append(primer_name)
        primer_dict[primer_name] = {}

        for assembly_file in os.listdir(assembly_dir_multi):
            assembly_path = os.path.join(assembly_dir_multi, assembly_file)
            org_name = assembly_file[:assembly_file.index(".")]
            amplicon_list.extend(ispcr(primer_path, assembly_path, MAX_AMPLICON_SIZE).split("\n"))
            #skip over the header lines and add sequences into a set to remove duplicate amplicons
            org_dict[org_name] = list({ampl for ampl in amplicon_list[1::2]})
            amplicon_list.clear() #clear list before next iteration

        #for each primer, add primer name as key and value as dict containing organism and resp. amplicons
        primer_dict[primer_name] = org_dict.copy()
        org_dict.clear() #clear dict before next iteration

    #write amplicons for each primer to multifasta file
    for primer, ampl_dict in primer_dict.items():
        with open(f"{primer}_amplicons.fna", "a") as fout:
            for org, ampl in ampl_dict.items():
                fout.write(f">{org}\n{ampl[0]}\n")

    #dictionary to keep track of org and concatenate the amplicons from different primer
    concat_ampl_dict = {}

    #perform MSA and concatenate the aligned amplicons   
    for primer in primer_list:        
        alined_file = msa(f"{primer}_amplicons.fna", f"{primer}_aln_ampl.efa")
        with open(alined_file, "r") as fin:
            line_list = fin.readlines()
            for line in line_list:
                if line.startswith(">"):
                    org = line.strip()
                    #initialize list for the first iteration
                    concat_ampl_dict[org] = [] if org not in concat_ampl_dict else concat_ampl_dict[org]
                    continue
                concat_ampl_dict[org].append(line.strip()) #aligned sequences

    #write the contents of dictionary into file
    with open("concat_aligned_ampl.efa", "a") as fout:
        for key, lis in concat_ampl_dict.items():
            fout.write(f"{key}\n")
            fout.writelines(lis)
            fout.write("\n")

    #generate newick string
    newick_str = newick_tree("concat_aligned_ampl.efa", "concat_newick_tree.txt")
    
    return newick_str

def msa(amplicon_file: str, out_file: str) -> str:
    muscle_out = subprocess.run([f"muscle -align {amplicon_file} -output {out_file}"], \
                   capture_output=True, \
                   text=True, \
                   shell=True, \
                   executable="/bin/bash"
                  )
    return out_file

def newick_tree(aligned_file: str, out_file:str) -> str:
    newick_str = subprocess.run([f"FastTree -nt {aligned_file}"],
                    capture_output=True, \
                    text=True, \
                    shell=True, \
                    executable="/bin/bash"
                   )
    return newick_str.stdout

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:
    """
    main function for calling functions to perform isPCR in three steps:
    1. Identify locations where primers would anneal to the target sequence
    2. Identify pairs of locations where two primers anneal close enough together and in the correct orientation for amplification to occur
    3. Extract the amplified sequence

    Args:
        primer_file: path to the pimer file
        assembly_file: path to the assembly file
        max_amplicon_size: maximum amplicon required

    Returns:
        amplicons that gets amplified in isPCR
    """
    sorted_good_hits = step_one(primer_file=primer_file, assembly_file=assembly_file)
    paired_hits = step_two(sorted_hits=sorted_good_hits, max_amplicon_size=max_amplicon_size)
    amplicons = step_three(hit_pairs=paired_hits, assembly_file=assembly_file)
    
    return amplicons

def step_one(primer_file: str, assembly_file: str) -> list[list[str]]:
    """
    to identify locations where primers would anneal to the target sequence

    Args:
        primer_file: file path of the primer file
        assembly_file: file path the assembly file

    Returns:
        list containing the filtered blast outputs
    """
    blast_output = call_blast(primer_file, assembly_file)
    filtered_blast_output = filter_blast(blast_output)
    
    return filtered_blast_output

def step_two(sorted_hits: list[list[str]], max_amplicon_size: int) -> list[tuple[list[str]]]:
    """
    identify pairs of locations where two primers anneal close enough together and in the correct 
    orientation for amplification to occur

    Args:
        sorted_hits: list of sorted and filtered blast outputs
        max_amplicon_size: maximum desired amplicon size

    Returns:
        a list of tuples with blast hit pairs which satisfy the condition of orientation and size
    """
    paired_hits = []
    paired_hits = find_amplicon_pairs(sorted_hits, max_amplicon_size, paired_hits)
    
    return paired_hits

def step_three(hit_pairs: list[tuple[list[str]]], assembly_file: str) -> str:
    """
    extracting amplicon sequences

    Args:
        hit_pairs: list of tuples with blast hit pairs which satisfy the condition of orientation and size
        assembly_file: file path for assembly file

    Returns:
        string with extracted amplicon sequences
    """
    bed_list = create_bed_file(hit_pairs)
    #create a string containing BED contents
    bed_content = ('\n').join(bed_list)

    amplicon = call_seqtk(assembly_file, bed_content)
    
    amplicon_list = amplicon.split("\n")
    
    for idx in range(len(bed_list)):
        if bed_list[idx].endswith("-"):
            #if amplicon in opposite orientation, reverse complement it
            amplicon_list[2*idx+1] = reverse_complement(amplicon_list[2*idx+1]) #skip over header lines

    return "\n".join(amplicon_list[:-1]) #skip the newline at the end

def call_blast(primer_file: str, assembly_file: str) -> str:
    """
    to blast primer sequence against assembly file to find sequence matches

    Args:
        primer_file: file path of the primer file
        assembly_file: file path the assembly file

    Returns:
        output of blast
    """
    blast_output = subprocess.run(["blastn", "-task", "blastn-short", "-query", primer_file, "-subject", assembly_file, "-outfmt", '6 std qlen'], \
                   capture_output=True, \
                   text=True
                   )
    return blast_output.stdout

def filter_blast(blast_output:str) -> list[list[str]]:
    """
    filter blast hits to extracts hits which match atleast a predefined threshold PCT_MATCH and store it as a list of strings

    Args:
        blast_output: output of blast

    Returns:
        list of sorted list of blast hits with match percent above PCT_MATCH
    """
    #$3 - percent macth identity
    #$4 - match length
    #$13 - query length
    filtered_blast_output = subprocess.run("awk '{if ($3 >= PCT_MATCH && $4 > 0.8*$13) print $0;}' | sort -k 9,10n", \
                                           capture_output=True, \
                                           text=True, \
                                           shell=True, \
                                           input=blast_output
                                           )
    blast_output_list = filtered_blast_output.stdout.split('\n')
    blast_output_fields = [i.split('\t') for i in blast_output_list[:-1]] #exclude last entry to get rid of unwanted newline

    return blast_output_fields

def find_amplicon_pairs(sorted_hits: list[list[str]], max_amplicon_size: int, paired_hits: list) -> list[tuple[list[str]]]:
    """
    Loop through all the sorted hits to check if any pair of hits satisfy conditions to make an amplicon
    1. Both primers anneal pointing towards one another
    2. Primers are sufficiently close to each other, set by max_amplicon_size

    Args:
        sorted_hits: list of sorted and filtered blast outputs
        max_amplicon_size (int): maximum desired amplicon size
        paired_hits (list): list to store hit pairs

    Returns:
        list of tuples with blast hit pairs which satisfy the condition of orientation and size
    """
    #primer[1] - query seqID
    #primer[0] - subject seqID
    #primer[8] and [9] - matched start and end position
    for primer1 in sorted_hits:
        for primer2 in sorted_hits:
            valid_amplicon_pair = ()            
            # We can skip comparing the same hits as they cannot make an amplicon
            if primer1 == primer2:
                continue
            else:
                # Check if both hits have the same sequence ID and primer IDs are not the same
                if primer1[1] == primer2[1] and primer1[0] != primer2[0]:
                    # Compare the 3' end of both primers and if their difference is less than amplicon size, they are valid pairs
                    if int(primer1[9]) < int(primer2[9]) and int(primer2[9]) - int(primer1[9]) < max_amplicon_size:
                        valid_amplicon_pair = (primer1, primer2)
                        paired_hits.append(valid_amplicon_pair)

    return paired_hits

def create_bed_file(hit_pairs: list[tuple[list[str]]]) -> list[str]:
    """
    create bed file using the filtered amplicon list

    Args:
        hit_pairs: list containing hit pairs

    Returns:
        BED content as string
    """
    bed_list = []
    #extract seqID, amplicon co-ordinate which is the 3' end position of the two primers
    for amplicon_pair in hit_pairs:
        primer1, primer2 = amplicon_pair
        if primer1[0].find("F") != -1:
            bed_list.append(f"{primer1[1]}\t{primer1[9]}\t{int(primer2[9])-1}\t+")
        else:
            bed_list.append(f"{primer1[1]}\t{primer1[9]}\t{int(primer2[9])-1}\t-")
    
    return bed_list

def call_seqtk(assembly_file: str, bed_content: str) -> str:
    """
    execute seqtk to extract sequence from assembly file using coordinates in bed file

    Args:
        assembly_file: path to assembly file
        bed_file: BED content

    Returns:
        string of amplicon sequences extracted from hit positions
    """
    amplicon_seq = subprocess.run(f'seqtk subseq -s {assembly_file} <(echo "{bed_content}" ; data/Vibrio_cholerae_N16961.bed | xargs)', \
                   shell=True, \
                   capture_output=True, \
                   text=True, \
                   executable="/bin/bash"
                   )
    return amplicon_seq.stdout

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

if __name__ == "__main__":
    tree = single_locus(primer_path, assembly_dir)
    print(tree)
    multi_tree = multi_loci(primer_dir_multi, assembly_dir_multi)
    print(multi_tree)

    