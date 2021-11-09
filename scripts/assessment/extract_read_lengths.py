import gzip
import os
import argparse




def count_fastqgz_reads_bases(fastq_gz):
    reads = 0
    bases = 0
    lengths = []
    with gzip.open(fastq_gz, 'r') as f:
        file_content = f.read().strip(b'\n').split(sep=b'\n')
        for i in range(0,len(file_content),4):
            if len(file_content) - i > 3:
                reads += 1
                bases += len(file_content[i+1])
                lengths.append(len(file_content[i+1]))
            else:
                print("WARNING: Truncated file: " + fastq_gz + " -> " + str(len(file_content)) + " " + str(i))

    return((lengths))

def count_fasta_length(fasta):
    fasta_file = open(fasta)
    header = fasta_file.readline()
    length = 0
    for line in fasta_file:
        length += len(line.strip('\n'))
    return(length)

def sample_sheet_to_used_ids_counts(sample_sheet):
    sh_file = open(sample_sheet)
    header = sh_file.readline()
    ids = []
    for line in sh_file:
        sp = line.strip().split("\t")
        count = count_fasta_length(sp[2])
        ids.append((sp[0].rjust(2,'0'),count))
    return(ids)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Track read destinations throughout the code')
    parser.add_argument('--sample_sheet', help='The sample sheet to load samples from',required=True)
    parser.add_argument('--results_dir', help='The base results directory',required=True)
    parser.add_argument('--output', help='an output file containing counts/coverage for samples at various stages of processing',required=True)

    args = parser.parse_args()
    sample_ids = sample_sheet_to_used_ids_counts(args.sample_sheet)
    print("samples: " + str(len(sample_ids)))

    out = open(args.output,"w")
    out.write("sample\treadLength\n")


    for i,count in sample_ids:

        try:
            #print(count_fastqgz_reads_bases(args.results_dir + "/length_filter/" + i + "_length_filtered.fq.gz"))
            for lengths in count_fastqgz_reads_bases(args.results_dir + "/filter_reads/" + i + "_filtered.fq.gz"):
                out.write(i + "\t" + str(lengths) + "\n")
            
            #filtered_coverage = count_fastqgz_reads_bases(args.results_dir + "/length_filter/" + i + "_length_filtered.fq.gz")
            #out.write(i + "\t" + "\t".join([str(x) for x in total_reads]) + "\t" + str(count) + "\t" + "\t".join([str(x) for x in guppy_coverage]) +
            #"\t" + "\t".join([str(x) for x in filtered_coverage]) + "\t" + "\t".join([str(x) for x in corrected_coverage]) + "\n")

        except:
            print("Couldn't read files for " + str(i))
