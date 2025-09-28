#!/usr/bin/env python3

import argparse
import pysam
from collections import defaultdict
from tqdm import tqdm

"""
Authors:
    Zoé Pochon: zoe.pochon@gmail.com

Description:
    This script processes a BAM file to remove C->T and G->A transitions, which are common damage patterns in ancient DNA,
    while retaining transitions at positions likely representing true mutations based on specific thresholds.

    The script uses two adjustable thresholds:
    - Transition frequency threshold (default 0.9): Positions with at least 90% transition frequency are considered true mutations.
    - Minimum depth threshold (default 3): Only positions with a minimum read depth of 3 are considered for mutation analysis, the other ones have transition removed.

    Usage:
        python maskDamage.py -b input.bam -o output.bam -r reference.fasta -t 0.9 -d 3
"""

# Parse arguments
parser = argparse.ArgumentParser(description="Process BAM file to remove C->T and G->A transitions unless they appear as true mutations.")
parser.add_argument("-b", "--bamfile", required=True, help="Path to the input BAM file.")
parser.add_argument("-o", "--output", required=True, help="Path to the output BAM file.")
parser.add_argument("-r", "--reference", required=True, help="Path to the reference FASTA file.")
parser.add_argument("-t", "--threshold", type=float, default=0.9, help="Threshold for transition percentage to consider as a true mutation.")
parser.add_argument("-d", "--min_depth", type=int, default=3, help="Minimum read depth required to consider a position a true mutation.")

args = parser.parse_args()

# Load the reference FASTA file
reference = pysam.FastaFile(args.reference)

# Initialize a dictionary to store counts for each position
transition_counts = defaultdict(lambda: {"total": 0, "transitions": 0})

# Counters for summary statistics
true_mutation_count = 0
masked_damage_count = 0

# Function to process each read and collect transition statistics
def count_transitions(bamfile):
    with pysam.AlignmentFile(bamfile, "rb") as bam:
        for read in tqdm(bam.fetch(), desc="Counting transitions"):
            if read.is_unmapped:
                continue

            ref_name = read.reference_name

            for query_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
                if ref_pos is None or query_pos is None:
                    continue  # skip insertions or clipping

                ref_base = ref_base.upper() if ref_base else reference.fetch(ref_name, ref_pos, ref_pos + 1).upper()
                read_base = read.query_sequence[query_pos]
                current_position = f"{ref_name}:{ref_pos}"

                transition_counts[current_position]["total"] += 1

                if (ref_base == "C" and read_base == "T") or (ref_base == "G" and read_base == "A"):
                    transition_counts[current_position]["transitions"] += 1

# Identify positions where transitions meet the threshold and minimum depth
def identify_true_mutations(threshold, min_depth):
    true_mutations = set()
    for position, counts in transition_counts.items():
        if counts["total"] >= min_depth and counts["transitions"] / counts["total"] >= threshold:
            true_mutations.add(position)
    return true_mutations

# Open input BAM and output BAM with header preserved
with pysam.AlignmentFile(args.bamfile, "rb") as bam, pysam.AlignmentFile(args.output, "wb", header=bam.header) as out_bam:
    count_transitions(args.bamfile)
    true_mutations = identify_true_mutations(args.threshold, args.min_depth)
    true_mutation_count = len(true_mutations)
    print(f"\nTransition analysis complete.")
    print(f"→ Identified {true_mutation_count} true mutation positions with depth ≥ {args.min_depth}.")

    for read in tqdm(bam.fetch(), desc="Processing reads"):
        if read.is_unmapped:
            continue

        read_seq = read.query_sequence
        read_qualities = read.query_qualities
        ref_pos = read.reference_start
        ref_name = read.reference_name
        updated_seq = list(read_seq)

        for query_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):
            if ref_pos is None or query_pos is None:
                continue

            current_position = f"{ref_name}:{ref_pos}"

            if current_position in true_mutations:
                continue

            ref_base = ref_base.upper() if ref_base else reference.fetch(ref_name, ref_pos, ref_pos + 1).upper()
            read_base = read.query_sequence[query_pos]

            if (ref_base == "C" and read_base == "T") or (ref_base == "G" and read_base == "A"):
                updated_seq[query_pos] = "N"
                masked_damage_count += 1

        read.query_sequence = "".join(updated_seq)
        read.query_qualities = read_qualities
        out_bam.write(read)

print()
print("----- Summary -----")
print(f"True mutation positions retained: {true_mutation_count}")
print(f"Likely damaged transitions masked with 'N': {masked_damage_count}")
print("-------------------")

