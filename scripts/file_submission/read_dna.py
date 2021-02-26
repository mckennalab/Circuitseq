#used in the convert_dna_to_fasta.sh, requires BioPython > 1.74, outputs an ugly fasta, fixed next
import sys
from Bio import SeqIO

count = SeqIO.convert(sys.argv[1], "snapgene", sys.argv[2], "fasta")
print("Converted %i records" % count)

