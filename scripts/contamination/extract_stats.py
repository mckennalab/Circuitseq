import pysam
import argparse


parser = argparse.ArgumentParser(description='Extract alignment statistics')
parser.add_argument('--bamfile', help='the bam file to load',required=True)
parser.add_argument('--output', help='our output stats file',required=True)
parser.add_argument('--sample', help='our sample',required=True)

args = parser.parse_args()

def bam_to_alignment_stats(bam_file):
    print(bam_file)
    samfile = pysam.AlignmentFile(bam_file, "rb")
    alignment_outcome_mapping = {}
    for i in range(0,10):
        alignment_outcome_mapping[i] = 0
    cnt = 0
    valid = 0
    total_bases = 0
    for read in samfile.fetch():
        cigar = read.cigartuples
        if cigar != None:
            valid += 1
            for tp in cigar:
                total_bases += tp[1]
                alignment_outcome_mapping[tp[0]] = alignment_outcome_mapping.get(tp[0],0) + tp[1]
        cnt += 1
    return((alignment_outcome_mapping,cnt,valid,total_bases))

output = open(args.output,"w")
output.write("index\tcount\tvalid_count\ttotal_bases\tM\tI\tD\tN\tS\tH\tP\tE\tX\tB\n")

stats = bam_to_alignment_stats(args.bamfile)
output.write(args.sample + "\t" + str(stats[1]) + "\t" + str(stats[2]) + "\t" + str(stats[3]) + "\t" + "\t".join([str(stats[0].get(x,0)) for x in range(0,10)]) + "\n")
output.close()
    