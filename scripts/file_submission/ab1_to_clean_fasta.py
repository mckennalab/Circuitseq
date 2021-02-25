#run as python ab1_to_fastq.py input output
from Bio import SeqIO
import sys

records = SeqIO.convert(sys.argv[1], "abi", sys.argv[2], "fastq")

with open(sys.argv[2],'r') as f:
    imported = []
    for i in f:
        imported.append(i.rstrip('\n'))
    sequence = imported[1]
f.close()

def findBiggestGap(s, ch):
    locs = ([i for i, letter in enumerate(s) if letter == ch])

    start = 0
    end = len(s)

    gap = 0

    for i in range(len(locs)-1):
        if locs[i+1] - locs[i] > gap:
            gap = locs[i+1] - locs[i]
            start = locs[i]
            end = locs[i+1]
    return(s[start+1:end])

largest_gap = findBiggestGap(sequence, 'N')

name  = (sys.argv[2].split('.'))[0]+'_clean.fasta'

out = open(name, 'w')
out.write('>'+(sys.argv[2].split('.'))[0]+'_clean\n')
out.write(largest_gap)
out.write('\n')
out.close()
