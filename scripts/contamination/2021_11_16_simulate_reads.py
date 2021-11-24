import pysam
import argparse
from shutil import copyfile
import gzip
import shutil
import os
import subprocess

# Create mixtures of contaminated reads from real, sequenced reads in rows:
# 1) another plasmid (plasmid 2 into plasmid 1)
# 2) the whole run, excluding the current well
# 3) the reverse of #1, plasmid 1 into plasmid 2
# 4) #2 but for plasmid 2

# in rows 5-8 we do the same BUT from purely simulated reads 

#parser = argparse.ArgumentParser(description='Create mixtures of contaminated reads')
#parser.add_argument('--bamfile', help='the bam file to load',required=True)
#parser.add_argument('--output', help='our output stats file',required=True)
#parser.add_argument('--sample', help='our sample',required=True)

# args = parser.parse_args()

def simulate_reads_old(plasmid, fastq_file, output_dir, reference_file):
    if not os.path.exists(output_dir):
        print("path doesn't exist. trying to make")
        os.makedirs(output_dir)
        
    # example commands:
    # read_analysis.py genome -i fastq_runid_97418955be55e4fdada56f8ae5845e80a66bf1b3_0.fastq -rg ../../results/rotated/49_rotated.fasta -o 49_sim -t 12
    # simulator.py genome -rg ../../results/rotated/49_rotated.fasta -dna_type circular -c 49_sim -n 5000
    # steps:
    # 1) cd into directory
    # 2) copy fastq.gz reads file
    # 3) unzip fastq file
    # 4) run analysis
    # 5) run simulator
    # 6) cleanup original reads file

    # 2 & 3
    fastq_new_file = output_dir + "/" + plasmid + ".fq"
    plasmid_name = "plasmid_sim_" + plasmid
    
    with gzip.open(fastq_file, 'rb') as f_in:
        with open(fastq_new_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    
    # 4
    output = subprocess.run(["read_analysis.py", "genome", "-i", fastq_new_file, "-rg", reference_file, "-o", plasmid_name], capture_output=True)
    print(output)
    
    # 5
    output = subprocess.run(["simulator.py", "genome", "-rg", reference_file,"-dna_type circular -n 5000" "-c", plasmid_name], capture_output=True)
    print(output)
    
    
def simulate_reads(plasmid, output_dir, reference_file):
    if not os.path.exists(output_dir):
        print("path doesn't exist. trying to make")
        os.makedirs(output_dir)
    
    # copy the reference and add the circular tag
    output_reference_file = output_dir + "/" + plasmid + ".fa"
    plasmid_length = 0
    
    with open(reference_file, 'r') as f_in:
        output_reference = open(output_reference_file,"w")
        header = f_in.readline().strip()
        output_reference.write(header + " circular=true\n")
        for line in f_in:
            output_reference.write(line)
            plasmid_length += len(line.strip())
        output_reference.close()
       
    cmd = 'python -m badread simulate --reference ' + output_reference_file + ' --cap_at_length --quantity 20M --length ' + str(plasmid_length) + ',5000 > ' + output_dir + '/reads.fastq' # | gzip > 
    print(cmd)
    mycmd=subprocess.getoutput(cmd)
    subprocess.getoutput('gzip ' + output_dir + '/reads.fastq')
    
for i in range(25,85):
    rjusted = str(i).rjust(2, '0')
    simulate_reads(str(i),
               "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/" + str(i) + "_test/",
               "/analysis/2021_08_26_PlasmidSeq_paper/results/rotated/" + rjusted + "_rotated.fasta")
for i in range(86,98):
    rjusted = str(i).rjust(2, '0')
    simulate_reads(str(i),
               "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/" + str(i) + "_test/",
               "/analysis/2021_08_26_PlasmidSeq_paper/results/rotated/" + rjusted + "_rotated.fasta")