
# ========================================================================================
#                          suffix array deduplicator
# ========================================================================================
#  attempt to remove duplicated segments that are 'too large' from the assembly
# ----------------------------------------------------------------------------------------
import numpy as np
from pydivsufsort import divsufsort, kasai
import argparse
import sys

parser = argparse.ArgumentParser(
    description="Attempt to remove duplicated regions from a plasmid"
)
parser.add_argument(
    "--plasmid_fasta",
    help="the input fasta containing an assembled reference",
    required=True,
)


parser.add_argument("--output_fasta", required=True, help="the output fasta file")

args = parser.parse_args()

ref = ""
reference_input = open(args.plasmid_fasta)
header = reference_input.readline()
for line in reference_input:
    ref += line.strip()
print(ref[0:20])


class CoverageMap:
    def __init__(
        self, string, strsuffix_array, string_lcp_array, minimum_merge_length=600
    ):
        self.string = string
        self.strsuffix_array = strsuffix_array
        self.string_lcp_array = string_lcp_array

        self.longest_coverage_segments = np.zeros(len(string), dtype=np.uintc)
        self.string_lcp_array_sorted_indexes = np.argsort(self.string_lcp_array)
        self.minimum_length_segments = []
        self.minimum_merge_length = minimum_merge_length

        for index, segment in enumerate(self.string_lcp_array_sorted_indexes):
            if segment + 1 < len(string):
                position1 = self.strsuffix_array[segment]
                position2 = self.strsuffix_array[segment + 1]
                seg_length = self.string_lcp_array[segment]

                self.longest_coverage_segments[
                    position1 : position1 + seg_length
                ] = index
                self.longest_coverage_segments[
                    position2 : position2 + seg_length
                ] = index

                if seg_length > minimum_merge_length:
                    self.minimum_length_segments.append(index)
        self.do_search = len(self.minimum_length_segments) > 0

    def find_co_segment_positions(self, segment_number):
        largest_segment_positions = np.where(
            self.longest_coverage_segments == segment_number
        )
        largest_segment_max_position = int(np.max(largest_segment_positions))
        largest_segment_max_start_position = int(
            np.max(largest_segment_positions) - (largest_segment_positions[0].size / 2)
        )
        largest_segment_min_position = int(np.min(largest_segment_positions))
        largest_segment_min_end_position = int(
            np.min(largest_segment_positions) + (largest_segment_positions[0].size / 2)
        )

        return (
            largest_segment_min_position,
            largest_segment_min_end_position,
            largest_segment_max_start_position,
            largest_segment_max_position,
            largest_segment_positions[0].size,
        )

    def grow_longest_segments(self, allowed_base_skips):
        assert self.minimum_merge_length > allowed_base_skips

        # find the longest segment
        counts = np.bincount(self.longest_coverage_segments)
        largest_segment = np.argmax(counts)
        absorbed_segments = []

        searching_for_extensions = True
        while searching_for_extensions:
            largest_segment_positions = self.find_co_segment_positions(largest_segment)
            searching_for_extensions = False
            # print(largest_segment_positions)
            # print(largest_segment_positions[2])
            left_slice_1 = self.longest_coverage_segments[
                max(
                    0, largest_segment_positions[0] - allowed_base_skips
                ) : largest_segment_positions[0]
            ]
            left_slice_2 = self.longest_coverage_segments[
                max(
                    0, largest_segment_positions[2] - allowed_base_skips
                ) : largest_segment_positions[2]
            ]
            left_intersect = np.intersect1d(left_slice_1, left_slice_2)
            # print("left intersect")
            # print(left_intersect)
            for left_overlap in left_intersect:
                if left_overlap in self.minimum_length_segments and not (
                    left_overlap in absorbed_segments
                ):
                    searching_for_extensions = True
                    absorbed_segments.append(left_overlap)
                    self.longest_coverage_segments[
                        max(
                            0, largest_segment_positions[0] - allowed_base_skips
                        ) : largest_segment_positions[0]
                    ] = largest_segment
                    self.longest_coverage_segments[
                        max(
                            0, largest_segment_positions[2] - allowed_base_skips
                        ) : largest_segment_positions[2]
                    ] = largest_segment
                    new_cover = np.where(self.longest_coverage_segments == left_overlap)
                    # print("COVER")
                    # print(len(new_cover))
                    self.longest_coverage_segments[new_cover] = largest_segment

            right_slice_1 = self.longest_coverage_segments[
                largest_segment_positions[1] : min(
                    largest_segment_positions[1] + allowed_base_skips, len(self.string)
                )
            ]
            right_slice_2 = self.longest_coverage_segments[
                largest_segment_positions[3] : min(
                    largest_segment_positions[3] + allowed_base_skips, len(self.string)
                )
            ]
            right_intersect = np.intersect1d(right_slice_1, right_slice_2)
            # print("right_intersect intersect")
            # print(right_intersect)
            for right_overlap in right_intersect:
                if right_overlap in self.minimum_length_segments and not (
                    right_overlap in absorbed_segments
                ):
                    searching_for_extensions = True
                    absorbed_segments.append(right_overlap)
                    self.longest_coverage_segments[
                        largest_segment_positions[1] : min(
                            largest_segment_positions[1] + allowed_base_skips,
                            len(self.string),
                        )
                    ] = largest_segment
                    self.longest_coverage_segments[
                        largest_segment_positions[3] : min(
                            largest_segment_positions[3] + allowed_base_skips,
                            len(self.string),
                        )
                    ] = largest_segment
                    new_cover = np.where(
                        self.longest_coverage_segments == right_overlap
                    )
                    # print("COVER")
                    # print(len(new_cover))
                    self.longest_coverage_segments[new_cover] = largest_segment

        self.largest_segment = largest_segment

    def subset_reference_to_repeat_free(self, allowed_base_skips):
        if self.do_search:
            self.grow_longest_segments(allowed_base_skips)
            segment_positions = self.find_co_segment_positions(self.largest_segment)
            segment1 = self.string[segment_positions[0] : segment_positions[2]]
            segment2 = self.string[segment_positions[3] : len(self.string)]
            segment3 = self.string[0 : segment_positions[0]]
            # print("DSFSFSDF")
            # print(segment1)
            # print(segment2)
            # print(segment3)
            print((segment1 + segment2 + segment3, segment_positions))
            return ((segment1 + segment2 + segment3, segment_positions))
        else:
            print((self.string, (0, 0, 0, 0, 0)))
            return((self.string, (0, 0, 0, 0, 0)))


def find_overlapping_fragments(input_string, minimum_length, spacing):
    lcps = []
    lcps_indexes = []
    best_ref = None
    best_rep_length = -1

    for i in range(0, len(input_string), spacing):
        # print("******* LOOP ******* ")
        rotation_string = input_string[i : len(input_string)] + input_string[0:i]
        # print(rotation_string)
        strsuffix_array = divsufsort(rotation_string)
        string_lcp_array = kasai(rotation_string, strsuffix_array)

        # sometimes small mismatches break up our fragments; try to extend over them
        cmap = CoverageMap(
            rotation_string, strsuffix_array, string_lcp_array, minimum_length
        )
        new_ref = cmap.subset_reference_to_repeat_free(300)

        print(new_ref)
        if new_ref[1][4] > best_rep_length:
            best_ref = new_ref
            best_rep_length = new_ref[1][4]

    return best_ref


lcps_plasmid = find_overlapping_fragments(ref, 500, 200)

output = open(args.output_fasta, "w")
output.write(
    header.strip() + "_" + "_".join([str(x) for x in lcps_plasmid[1:6]]) + "\n"
)
output.write(lcps_plasmid[0] + "\n")
