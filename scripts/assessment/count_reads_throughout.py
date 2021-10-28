import gzip
import os
import argparse

def count_fastq_reads_bases(fastq):
    reads = 0
    bases = 0
    with open(fastq) as f:
        file_content = f.read().split()
        for i in range(0,len(file_content),4):
            reads += 1
            bases += len(file_content[i+1])

    return((reads,bases))

def count_fastq_gz_dir(directory):
    reads = 0
    bases = 0
    for fl in os.listdir(directory):
        if ".fastq.gz" in fl:
            try:
                (read,base) = count_fastqgz_reads_bases(directory + "/" + fl)
            except:
                print("failed on file " + fl)
            reads += read
            bases += base
    return((reads,bases))


def count_fasta_length(fasta):
    fasta_file = open(fasta)
    header = fasta_file.readline()
    length = 0
    for line in fasta_file:
        length += len(line.strip('\n'))
    return(length)

def count_fastagz_reads_bases(fastq_gz):
    reads = 0
    bases = 0
    with gzip.open(fastq_gz, 'r') as f:
        file_content = f.read().strip(b'\n').split(sep=b'\n')
        for line in file_content:
            if line.startswith(b'>'):
                reads += 1
            else:
                bases += len(line)

    return((reads,bases))

def count_fastqgz_reads_bases(fastq_gz):
    reads = 0
    bases = 0
    with gzip.open(fastq_gz, 'r') as f:
        file_content = f.read().strip(b'\n').split(sep=b'\n')
        for i in range(0,len(file_content),4):
            if len(file_content) - i > 3:
                reads += 1
                bases += len(file_content[i+1])
            else:
                print("WARNING: Truncated file: " + fastq_gz + " -> " + str(len(file_content)) + " " + str(i))

    return((reads,bases))

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
    out.write("sample\ttotalReads\ttotalCounts\tlength\tguppyReads\tguppyBases\tfilteredReads\tfilteredBases\tcorrectedReads\tcorrectedBases\n")

    total_reads = count_fastq_gz_dir(args.results_dir + "/guppy/basecalling/pass/")

    for i,count in sample_ids:
        try:
            guppy_coverage = count_fastq_gz_dir(args.results_dir + "/guppy_demultiplex/saved_data/barcode" + i )
            filtered_coverage = count_fastqgz_reads_bases(args.results_dir + "/length_filter/" + i + "_length_filtered.fq.gz")
            corrected_coverage = count_fastagz_reads_bases(args.results_dir + "canu/" + i + "_canu_correct/reads.correctedReads.fasta.gz")
            out.write(i + "\t" + "\t".join([str(x) for x in total_reads]) + "\t" + str(count) + "\t" + "\t".join([str(x) for x in guppy_coverage]) +
            "\t" + "\t".join([str(x) for x in filtered_coverage]) + "\t" + "\t".join([str(x) for x in corrected_coverage]) + "\n")

        except:
            print("Couldn't read files for " + str(i))
