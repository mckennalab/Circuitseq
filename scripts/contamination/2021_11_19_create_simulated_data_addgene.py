import pysam
import argparse
import random 
import gzip

# mix reads to create synthetically contaminated samples
sample_sheet = open("computational_96_well_mixture_addgene.txt","w")
sample_sheet.write("position\tsample\treference\tsangers\treads\n")

# take the first four plasmids, use those as the intended plasmids
# then for odd numbered rows, mix with 96 - (plasmid number) at 
# ratios 0%, 1%, 2%, 3%, 5%, 10%, 15%, 20%, 25%, 30%, 40%, 50%
# for the even rows do the same, but draw reads randomly from all
# other wells on the plate.

# base_plasmids = [1,2,3,4]
replicates = ['a','b','c','d','e']

top_15_addgene = [12260,12259,48138, 52961,42230,62988,8454,12253,12251,1864,10878,21915,52963,48140,48139]
blue_flames = [71720,52535,51692,87360,14723,24595, 30313, 12301, 66801,68122,
               47327,10808,40342,32751,45161,14873, 107177,47504, 81070,12298,
               36939,35617,14436,80925,59946,110060,26722, 26646, 49792,66810,
               31485,32515,1654, 52213,43796,64073, 20960, 36325, 20737,11427,
               31815,14080,16337,47443,31367,52925, 83467, 45605, 27052,8381,
               65974,35175,57822,36412,22945,21232, 37120, 16542, 10953,70219,
               29435,47549,8449, 61463,27340,14129, 49172, 104991,60415,18110,
               11916,40729,75127,73501,1015, 11181, 50473, 39196, 26678,63890,11908]

standards = [31815,16337,12251,49792,52961]
# mixing_single_plasmids = [5,6,7,8]
mixing_props = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99]

# /analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/addgene_simulated_data/12251_test/reads.fastq.gz
base_output = "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/simulated_mixed/"
#simlulated_reads = "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/simulated_reads/"
simlulated_reads = "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/addgene_simulated_data/"
simlulated_reads_ending = "_test/reads.fastq.gz"

# load up the previous sample sheet to get references
known_refs = open("known_plasmids.tsv")
header = known_refs.readline()
sample_to_ref = {}
for line in known_refs:
    sp = line.strip().split("\t")
    sample_to_ref[int(sp[0])] = sp[2]

# helper functions
def count_reads(file):
    cnt = 0
    with pysam.FastxFile(file) as fh:
        for entry in fh:
            cnt += 1
    return(cnt)

def sample_reads(file,number_to_sample,reads=None):
    if reads == None:
        reads = all_reads(file)
    return(random.sample(reads,number_to_sample))

def all_reads(file):
    all_reads = []
    with pysam.FastxFile(file) as fh:
            for entry in fh:
                all_reads.append(entry)
    return(all_reads)

print("reading all background reads...")
background_all_reads = all_reads("all_other_addgene_reads.gz")

cnt = 0
for index, source_plasmid in enumerate(standards):
    source_name = source_plasmid # str(source_plasmid).rjust(2, '0')
    source_file = simlulated_reads + str(source_plasmid) + simlulated_reads_ending
    original_read_count = count_reads(source_file)
                    
    for replicate in replicates:  
        for mixing_rate in mixing_props:
            # create the many sample mixing data set
            sample_name = str(source_name) + "_common_at_" + str(mixing_rate).replace(".","p") + "_rep_" + replicate
            print(sample_name)
            cnt +=1 
            sample_sheet.write(str(cnt) + "\t" + sample_name + "\t" + simlulated_reads + str(source_plasmid) + "_test/" + str(source_plasmid)  + ".fa\t\tsimulated_mixed/" + sample_name + ".fastq.gz\n")
            # gzip.open('/home/joe/file.txt.gz', 'wb')
            new_output =  gzip.open(base_output + sample_name + ".fastq.gz","wb")

            sampled_reads_contam = sample_reads("dummy",int((mixing_rate * original_read_count)/2),background_all_reads)
            if mixing_rate > 0:
                sampled_reads_original = sample_reads(source_file,int((original_read_count * (1.0 - mixing_rate))/2))
            else:
                sampled_reads_original = all_reads(source_file)

            sampled_reads_contam.extend(sampled_reads_original)
            random.shuffle(sampled_reads_contam)

            for entry in sampled_reads_contam:
                new_output.write(("@" + entry.name + "\n" + entry.sequence + "\n+\n" + entry.quality + "\n").encode())
            new_output.close()
        
sample_sheet.close()            
            
            
        