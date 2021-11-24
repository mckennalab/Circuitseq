import pysam
import argparse
import random 
import gzip

# mix reads to create synthetically contaminated samples
sample_sheet = open("computational_96_well_mixture.txt","w")
sample_sheet.write("position\tsample\treference\tsangers\treads\n")

# take the first four plasmids, use those as the intended plasmids
# then for odd numbered rows, mix with 96 - (plasmid number) at 
# ratios 0%, 1%, 2%, 3%, 5%, 10%, 15%, 20%, 25%, 30%, 40%, 50%
# for the even rows do the same, but draw reads randomly from all
# other wells on the plate.

base_plasmids = [1,2,3,4]
replicates = ['a','b','c','d','e']
# mixing_single_plasmids = [5,6,7,8]
mixing_props = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99]

base_output = "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/simulated_mixed/"
simlulated_reads = "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/simulated_reads/"
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
background_all_reads = all_reads("all_reads_minus_1_2_3_4.fastq.gz")

cnt = 0
for index, source_plasmid in enumerate(base_plasmids):
    source_name = str(source_plasmid).rjust(2, '0')
    source_file = simlulated_reads + str(source_plasmid) + simlulated_reads_ending
    original_read_count = count_reads(source_file)
                    
    for replicate in replicates:  
        for mixing_rate in mixing_props:
            # create the many sample mixing data set
            sample_name = source_name + "_common_at_" + str(mixing_rate).replace(".","p") + "_rep_" + replicate
            print(sample_name)
            cnt +=1 
            sample_sheet.write(str(cnt) + "\t" + sample_name + "\t" + sample_to_ref[source_plasmid] + "\t\tsimulated_mixed/" + sample_name + ".fastq.gz\n")
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
            
            
        