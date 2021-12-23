
# ========================================================================================
#                          suffix array deduplicator
# ========================================================================================
#  attempt to remove duplicated segments that are 'too large' from the assembly
# ----------------------------------------------------------------------------------------
import itertools
import numpy as np
import argparse
import sys
import math
from pydivsufsort import divsufsort, kasai

parser = argparse.ArgumentParser(
    description="Attempt to remove duplicated regions from a plasmid"
)
parser.add_argument(
    "--plasmid_fasta",
    help="the input fasta containing an assembled reference",
    required=True,
)

parser.add_argument("--output_fasta", required=True, help="the output fasta file")

args = parser.parse_args()


def matrix(seq1, seq2, match_score=2, gap_cost=1,min_bandwidth=100):
    assert(len(seq1) >= len(seq2))
    
    H = np.zeros((len(seq1), len(seq2)), np.int)

    for i in range(0,H.shape[0]):
        for j in range(0,H.shape[1]):
            match = H[i - 1, j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else - match_score)
            delete = H[i - 1, j] - gap_cost
            insert = H[i, j - 1] - gap_cost
            H[i, j] = max(match, delete, insert, 0)
    return H

def traceback(matrix, a,b):
    max_index = np.unravel_index(matrix.argmax(), matrix.shape)
    
    index_i = max_index[0]
    index_j = max_index[1]
    
    substr_a = a[index_i]
    substr_b = b[index_j]
    
    while index_i > 0 and index_j > 0 and matrix[index_i][index_j] > 0:
        transitions = np.argmax([matrix[index_i-1,index_j],matrix[index_i-1,index_j-1],matrix[index_i,index_j-1]])
        if transitions == 0:
            substr_a += a[index_i - 1]
            substr_b += '-'
            index_i -= 1
        elif transitions == 1:
            substr_a += a[index_i - 1]
            substr_b += b[index_j - 1]
            index_i -= 1
            index_j -= 1
        elif transitions == 2:
            substr_a += '-'
            substr_b += b[index_j - 1]
            index_j -= 1
        else:
            raise NameError("did something dumb")
    return((matrix[max_index[0],max_index[1]],substr_a[::-1],substr_b[::-1],max_index,(index_i,index_j)))

def load_plasmid(fasta):
    pl_file = open(fasta)
    header = pl_file.readline()
    sequence = ""
    for line in pl_file:
        sequence += line.strip().upper()
    return((header,sequence))

def create_final_plasmid(sequence1,sequence2,subseq1,subseq2,score,minimum_length=1000,minumum_score_prop=1.8):
    # create a plasmid that's the first segment plus the second, which will preserve the alignment orientation
    full = sequence1 + sequence2
    if max(len(subseq1),len(subseq2)) > minimum_length and float(score)/max(len(subseq1),len(subseq2)) > minumum_score_prop:
        # we drop the smaller seqment from the plasmid
        no_gap_seg1 = "".join([x if x != '-' else "" for x in subseq1])
        no_gap_seg2 = "".join([x if x != '-' else "" for x in subseq2])
        no_gap_seg1_pos = full.index(no_gap_seg1)
        no_gap_seg2_pos = full.index(no_gap_seg2)
        
        if len(no_gap_seg1) > len(no_gap_seg2):
            print(len(full[0:no_gap_seg2_pos] + full[no_gap_seg2_pos+len(no_gap_seg2):len(full)]))
            return(full[0:no_gap_seg2_pos] + full[no_gap_seg2_pos+len(no_gap_seg2):len(full)])
        else:
            print(len(full[0:no_gap_seg1_pos] + full[no_gap_seg1_pos+len(no_gap_seg1):len(full)]))
            return(full[0:no_gap_seg1_pos] + full[no_gap_seg1_pos+len(no_gap_seg1):len(full)])
    else:
        print(len(full))
        return(full)
        
    

fasta_file = args.plasmid_fasta
fasta_header = load_plasmid(fasta_file)
fasta = fasta_header[1]
header = fasta_header[0]
output = open(args.output_fasta, "w")
output.write(header)

strsuffix_array = divsufsort(fasta)
string_lcp_array = kasai(fasta, strsuffix_array)

max_pos = np.argmax(string_lcp_array)
max_start = strsuffix_array[max_pos]
max_start2 = strsuffix_array[max_pos+1]

segment1 = fasta[max_start:max_start2]
segment2 = fasta[max_start2:len(fasta)] + fasta[0:max_start]

if max_start > max_start2:
    segment1 = fasta[max_start2:max_start]
    segment2 = fasta[max_start:len(fasta)] + fasta[0:max_start2]
    
if len(segment1) > len(segment2):
    mat = matrix(segment1,segment2)
    print("traceback")
    tb = traceback(mat,segment1,segment2)
    full = create_final_plasmid(segment1,segment2,tb[1],tb[2],tb[0])
    output.write(full + "\n")
else:
    mat = matrix(segment2,segment1)
    print("traceback")
    tb = traceback(mat,segment2,segment1)
    full = create_final_plasmid(segment2,segment1,tb[1],tb[2],tb[0])
    output.write(full + "\n")

output.close()
