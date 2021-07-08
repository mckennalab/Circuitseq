import numpy as np
from pydivsufsort import divsufsort, kasai
import argparse
import sys

parser = argparse.ArgumentParser(
    description="Attempt to align an assembled plasmid and a reference"
)
parser.add_argument(
    "--assembly",
    help="the input fasta containing an assembled reference",
    required=True,
)
parser.add_argument(
    "--reference",
    help="the input reference fasta containing an assembled reference",
    required=True,
)
parser.add_argument(
    "--reference_out",
    help="the reference hopefully aligned to the assembly",
    required=True,
)
parser.add_argument(
    "--assembly_out",
    help="the assembly hopefully aligned to the reference",
    required=True,
)
args = parser.parse_args()


def rev_base(base):
    if base == "A":
        return "T"
    if base == "C":
        return "G"
    if base == "G":
        return "C"
    if base == "T":
        return "A"
    return "N"


def reverse_comp(seq):
    return "".join([rev_base(b) for b in seq][::-1])


def find_longest_prefix(string1, string2):
    max_suffix = 0
    max_index = 0

    joined = string1 + "Z" + string2
    strsuffix_array = divsufsort(joined)
    string_lcp_array = kasai(joined, strsuffix_array)

    for i in range(0, len(strsuffix_array) - 1):
        is_ref1 = strsuffix_array[i] > len(ref)
        is_ref2 = strsuffix_array[i + 1] > len(ref)
        if string_lcp_array[i] > max_suffix and is_ref1 != is_ref2:
            max_suffix = string_lcp_array[i]
            max_index = i

    if strsuffix_array[max_index] < strsuffix_array[max_index + 1]:
        return (
            max_suffix,
            max_index,
            strsuffix_array[max_index],
            strsuffix_array[max_index + 1] - (1 + len(string1)),
        )
    return (
        max_suffix,
        max_index,
        strsuffix_array[max_index + 1],
        strsuffix_array[max_index] - (1 + len(string1)),
    )


ref = ""
reference_input = open(args.reference)
ref_header = reference_input.readline().strip()
for line in reference_input:
    ref += line.strip().upper()

assembly = ""
assembly_input = open(args.assembly)
assembly_header = assembly_input.readline().strip()
for line in assembly_input:
    assembly += line.strip().upper()


fwd = find_longest_prefix(ref, assembly)
rev_comp = find_longest_prefix(ref, reverse_comp(assembly))

ref_new = ""
assembled_new = ""

if fwd[0] > rev_comp[0]:
    print("FWD")
    ref_new = ref[fwd[2] : len(ref)] + ref[0 : fwd[2]]
    assembled_new = assembly[fwd[3] : len(assembly)] + assembly[0 : fwd[3]]
else:
    print("RC")
    ref_new = ref[rev_comp[2] : len(ref)] + ref[0 : rev_comp[2]]
    assembled_new = (
        reverse_comp(assembly)[rev_comp[3] : len(assembly)]
        + reverse_comp(assembly)[0 : rev_comp[3]]
    )

output_reference = open(args.reference_out, "w")
output_reference.write(ref_header + "\n")
output_reference.write(ref_new + "\n")
output_reference.close()

output_assembly = open(args.assembly_out, "w")
output_assembly.write(assembly_header + "\n")
output_assembly.write(assembled_new + "\n")
output_assembly.close()
