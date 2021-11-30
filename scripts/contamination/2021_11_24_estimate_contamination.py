import pysam
import argparse
import math
import numpy as np
import argparse
import cProfile
import random

parser = argparse.ArgumentParser(description='Estimate contamination rates for a plasmid')
parser.add_argument('--bamfile', help='the aligned reads bam file',required=True)
parser.add_argument('--reference', help='our reference file',required=True)
parser.add_argument('--sample', help='the sample name',required=True)
parser.add_argument('--plasmid_mix_in_rate', help='the estimated rate of plasmid mix-in',required=True,type=float)
parser.add_argument('--sample_estimate', help='the sample contamination, as a single line file',required=True)
parser.add_argument('--sample_range', help='the sample contamination estimations over the whole spectrum',required=True)

args = parser.parse_args()

def code_to_cigar_char(code):
    if code == 0:
        return("M")
    elif code == 1:
        return("I")
    elif code == 2:
        return("D")
    elif code == 3:
        return("N")
    elif code == 4:
        return("S")
    elif code == 5:
        return("H")
    elif code == 6:
        return("P")
    elif code == 7:
        return("=")
    elif code == 8:
        return("X")
    elif code == 9:
        return("B")

    
    
def aligned_reference_read(reference,read,start,cigar_tuples,quals,debug=False):
    '''
    Take the alignment and convert it into a simple aligned/unaligned call for
    every reference/read basepair within the read
    '''
    reference_position = start
    read_position = 0
    return_alignments = []
    return_quals = []
    #print(reference)
    #print(read)
    
    # tuples are [event_type,event_length]
    try:
        for cigar_tup in cigar_tuples:
            if cigar_tup[0] == 0: # match/mismatch
                for i in range(0,cigar_tup[1]):
                    if reference[reference_position] == read[read_position]:
                        return_alignments.append(True)
                    else:
                        return_alignments.append(False)
                    return_quals.append(quals[read_position])
                    reference_position += 1
                    read_position += 1


            elif cigar_tup[0] == 1: # insertion
                return_alignments.extend([False] * cigar_tup[1])
                return_quals.extend(quals[read_position:read_position + cigar_tup[1]])
                read_position += cigar_tup[1]

            elif cigar_tup[0] == 2: # deletion
                reference_position += cigar_tup[1] # not sure what we should do with for our contamination scores

            elif cigar_tup[0] == 3: # ref skip? 
                raise NameError("We dont support skips")

            elif cigar_tup[0] == 4: # soft clip
                for i in range(0,cigar_tup[1]):
                    return_alignments.append(False)
                    return_quals.append(quals[read_position])
                    # print("SOFT " + str(len(reference)) + "\t" + str(len(read)) + "\t" + str(reference_position) + "\t" + str(read_position) + "\t" + str(i)) 
                    read_position += 1

            elif cigar_tup[0] == 5: # hard clip
                for i in range(0,cigar_tup[1]):
                    return_alignments.append(False)
                    return_quals.append(chr(12))
                    # read_position += 1

            elif cigar_tup[0] == 6: # pad?
                raise NameError("We dont support pads")

            elif cigar_tup[0] == 7: # equal?
                raise NameError("We dont support equals")

            elif cigar_tup[0] == 8: # diff?
                raise NameError("We dont support diffs")

            elif cigar_tup[0] == 9: # back?
                raise NameError("We dont support backs")

            if debug:
                print(code_to_cigar_char(cigar_tup[0]) + "-" + str(cigar_tup[1]) + "\tref=" + str(len(reference)) + "\tread=" + str(len(read)) + "\trefpos=" + str(reference_position) + "\treadpos=" + str(read_position) + "\t" + str(i)) 


        return(return_alignments)
    except:
        print("failed on " + read + " ------------and------------- " + reference)

# aligned_reference_read("AAAAAAA","AAAT",0,[(0,2),(2,3),(0,2)])

our_match_inflation = args.plasmid_mix_in_rate
precomputed_error_rate_contamination_scores = {}
for i in range(0,100):
    for j in range(0,1000):
        for k in [True,False]:
            error_rate = math.pow(10,-1.0 * (i/10))
            contamination_score = j / 1000
            if j == 1000:
                contamination_score = 0.9999999999999999
            elif j == 0:
                contamination_score = 0.0000000000000001
            if k:
                precomputed_error_rate_contamination_scores[(error_rate,j/1000,k)] = math.log( ((1.0 - contamination_score) * (1.0 - error_rate)) + (contamination_score * (.25 )))
            else:
                precomputed_error_rate_contamination_scores[(error_rate,j/1000,k)] = math.log( ((1.0 - contamination_score) * error_rate) + (contamination_score * (0.75 )))
           
        
def match_score(bases,error_rates,matches,contamination_score_input,aligned):    
    return(np.sum([precomputed_error_rate_contamination_scores[(error_rates[i],contamination_score_input,matches[i])] for i in range(0,len(bases))]))


def update_scores(scores, score_bins, bases,quals,matches,norm):
    subsampled_positions = random.sample(range(len(bases)), 500)
    bases = [bases[i] for i in subsampled_positions]
    quals =[quals[i] for i in subsampled_positions]
    matches = [matches[i] for i in subsampled_positions]
    error_rates = [math.pow(10,-1.0 * (q/10)) for q in quals]
    for i in range(0,len(scores)):
        scores[i] += match_score(bases,error_rates,matches,score_bins[i],norm)
    
def bam_to_alignment_stats(reference, bam_file,bins=100):

    log_contamination_bins = np.zeros(bins) # [0 for x in range(0,bins)]
    contamination_representations = [x/bins for x in range(0,bins)]
    
    print(bam_file)
    
    samfile_pre = pysam.AlignmentFile(bam_file, "rb")
    lens = []
    for read in samfile_pre.fetch():
        
        sequence  = read.query_sequence
        qualities = read.query_qualities
        
        if sequence != None and not read.is_secondary and len(sequence) > 500:
            lens.append(len(sequence))
    
    max_length = max(lens)
    min_length = min(lens)
    print("min read length " + str(min_length) + " min read length " + str(max_length))# + " total " + str(len(lens))
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    total_reads = 0
    unprocessed = 0
    matches = 0
    total = 0
    good_reads = 0
    total_filter = 0
    for read in samfile.fetch():
        total_reads += 1
        if total_reads % 100 == 0:
            print("Processed = " + str(total_reads) + " unprocessed = " + str(unprocessed) + " max " + str(np.argmax(log_contamination_bins)))
            # break
        sequence  = read.query_sequence
        qualities = read.query_qualities
        
        if sequence == None or read.is_secondary or len(sequence) < 500:
            unprocessed += 1
        elif read.is_unmapped:
            aligned_set = [True if random.random() < (.25) else False for index, x in enumerate(read.query_sequence)]
            matches     += sum(aligned_set)
            total += len(aligned_set)
            total_filter += 1
            norm = math.log(max(0.0001,1.0 - (len(sequence) - min_length)/(max_length)))
            update_scores(log_contamination_bins, contamination_representations, sequence, qualities, aligned_set,norm)
        else:
            position    = read.get_reference_positions()[0]
            cigar       = read.cigartuples
            aligned_set = aligned_reference_read(reference,sequence,position,cigar,qualities)
            matches     += sum(aligned_set)
            total_filter += 1
            if sum(aligned_set) / len(sequence) > .7:
                good_reads += 1
            total += len(aligned_set)
            norm = math.log(max(0.0001,1.0 - (len(sequence) - min_length)/(max_length)))
            update_scores(log_contamination_bins, contamination_representations, sequence, qualities, aligned_set,norm)
    
    print(matches)
    print(total)
    print(good_reads)
    print(total_filter)
    return((log_contamination_bins,contamination_representations))
            
ref = ""
plasmid_ref = open(args.reference)
header = plasmid_ref.readline()
for line in plasmid_ref:
    ref += line.strip().upper()

bin_count = 100
# cProfile.run('bins = bam_to_alignment_stats(ref,args.bamfile,bin_count)')
bins = bam_to_alignment_stats(ref,args.bamfile,bin_count)

single_line = open(args.sample_estimate,"w")
single_line.write(args.sample + "\t" + str(np.argmax(bins[0])) + "\n")
single_line.close()


raw_probs = [math.exp(x - np.max(bins[0])) for x in bins[0]]
probs = [x/sum(raw_probs) for x in raw_probs]

all_lines = open(args.sample_range,"w")
for i in range(0,bin_count):
    all_lines.write(args.sample + "\t" + str(i/bin_count) + "\t" + str(probs[i]) + "\n")
all_lines.close()