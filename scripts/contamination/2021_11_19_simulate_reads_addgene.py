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
       
    cmd = 'python -m badread simulate --reference ' + output_reference_file + ' --cap_at_length --quantity 20M --emp_fragmentsemp_fragments /analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/distro_plasmid_lengths_normalized_2021_11_19.txt,' + output_reference_file + ' --length ' + str(plasmid_length) + ',5000 > ' + output_dir + '/reads.fastq' # | gzip > 
    print(cmd)
    mycmd=subprocess.getoutput(cmd)
    
    outstats = open(output_dir + "/read_sizes.txt","w")
    instats = open(output_dir + '/reads.fastq')
    for index,line in enumerate(instats):
        if index % 4 == 1:
            outstats.write(str(len(line)) + "\n")
    outstats.close()
    
    subprocess.getoutput('gzip ' + output_dir + '/reads.fastq')

top_15_addgene = [12260,12259,48138, 52961,42230,62988,8454,12253,12251,1864,10878,21915,52963,48140,48139]
blue_flames = [71720,52535,51692,87360,14723,24595, 30313, 12301, 66801,68122,
               47327,10808,40342,32751,45161,14873, 107177,47504, 81070,12298,
               36939,35617,14436,80925,59946,110060,26722, 26646, 49792,66810,
               31485,32515,1654, 52213,43796,64073, 20960, 36325, 20737,11427,
               31815,14080,16337,47443,31367,52925, 83467, 45605, 27052,8381,
               65974,35175,57822,36412,22945,21232, 37120, 16542, 10953,70219,
               29435,47549,8449, 61463,27340,14129, 49172, 104991,60415,18110,
               11916,40729,75127,73501,1015, 11181, 50473, 39196, 26678,63890,11908]

print(len(blue_flames))
all_plasmids = top_15_addgene + blue_flames

for index, i in enumerate(all_plasmids):
    rjusted = str(index).rjust(2, '0')
    simulate_reads(str(i),
               "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/addgene_simulated_data/" + str(i) + "_test/",
               "/analysis/2021_08_26_PlasmidSeq_paper/scripts/contamination/plasmids_blue_flame/plasmid_download_" + str(i) + ".fa")
