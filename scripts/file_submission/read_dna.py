import sys
from Bio import SeqIO

count = SeqIO.convert(sys.argv[1], "snapgene", sys.argv[2], "fasta")
print("Converted %i records" % count)

