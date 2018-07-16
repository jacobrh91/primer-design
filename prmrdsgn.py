#!/usr/bin/env python3

#  Written by: Hans Müller Paul and Jacob Heldenbrand

# IMPORT LIBRARIES

import datetime
import itertools
import argparse
import math
import sys
import statistics
import json
from os import path

# DEFINING VARIABLES

parser = argparse.ArgumentParser()

parser.add_argument('-i', dest='i', required=True, #metavar='', defining metavar is causing an error for some reason
                    help='path to input file in FASTA format'
                    )
parser.add_argument('-n', '--number', metavar='', dest='n', type=int, default=5,
                    help='number of candidate primer pairs to pick, default = 5'
                    )
parser.add_argument('-e', '--extension', metavar='', dest='e', type=int, default=100,
                    help='number of bases from the start and end of the sequence to look for primers in, default = 100'
                    )
parser.add_argument('-s', '--short', metavar='', dest='s', type=int, default=20,
                    help='shortest acceptable primer, default = 20'
                    )
parser.add_argument('-l', '--long', metavar='', dest='l', type=int, default=30,
                    help='longest acceptable primer, default = 30'
                    )
parser.add_argument('-m', '--mintemp', metavar='', dest='m', type=float, default=50,
                    help='min Tm in celsius, default = 55'
                    )
parser.add_argument('-x', '--maxtemp', metavar='', dest='x', type=float, default=65,
                    help='max Tm in celsius, default = 62'
                    )
parser.add_argument('-M', '--mingc', metavar='', dest='M', type=float, default=35, 
                    help='min GC percentage, default = 40'
                    )
parser.add_argument('-X', '--maxgc', metavar='', dest='X', type=float, default=65, 
                    help='max GC percentage, default = 60'
                    )
parser.add_argument('-D', '--tmdiff', metavar='', dest='D', type=float, default=0.5,
                    help='accepted TM difference to form primer pair'
                    )
parser.add_argument('-o', '--output', metavar='', dest='o', default='output.txt',
                    help='path to output file in FASTA format'
                    )
parser.add_argument('-v', '--verbose', action='store_true',
                    help='prints each step of each iteration (for debugging)'
                    )

args = parser.parse_args()


class Primer:
    def __init__(self, sequence):
        self.sequence = sequence
        self.GC_percentage = self.__calculate_GC()
        self.Tm = self.__calculate_Tm()

    def __calculate_GC(self):
        upper_seq = self.sequence.upper()
        gc_count = upper_seq.count('G') + upper_seq.count('C')
        gc_fraction = float(gc_count) / len(self.sequence)
        return 100 * gc_fraction

    def __calculate_Tm(self):
        N = len(self.sequence)
        melt_temp = 81.5 + 0.41 * self.GC_percentage - (675 / N)
        return melt_temp

    def get_primer_elements(self):
        return {"sequence": self.sequence,
                "GC": self.GC_percentage,
                "Tm": self.Tm
                }


def create_reverse_complement(input_sequence):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    bases = list(input_sequence)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    return bases


def get_primers(input_sequence):
    """
     Given an input sequence, return a list of all possible Primer objects
    """
    # This was len(input_sequence[:args.e]), which is equivalent to args.e
    length = args.e
    primer_sequence_list = [input_sequence[i:j + 1] for i in range(length) for j in range(i + args.s, i + args.l)]
    # Turn primer sequences into Primer objects
    primer_list = [Primer(primer) for primer in primer_sequence_list]
    return primer_list


def filter_primers(list_of_primers):
    filtered_list = list_of_primers.copy()
    for primer in list_of_primers:
        if primer.GC_percentage < args.M or primer.GC_percentage > args.X or primer.Tm < args.m or primer.Tm > args.x:
            filtered_list.remove(primer)
    return filtered_list


def main():

    if args.verbose:
        print(args)

    # IMPORTING FASTA FILE
    if path.exists(args.i):
        with open(args.i, 'r') as f:
            contents = f.read()
            seqlist = contents.split()
        name_tag = seqlist[0]
        seq = ''.join(seqlist[1:])
        if args.verbose:
            print(f'file name is {name_tag} and sequence is {seq}')
    else:
        sys.exit('file does not exist')

    # Lists of Primer objects
    forward_primers = get_primers(seq)
    reverse_primers = get_primers(create_reverse_complement(seq))

    filtered_forward_primers = filter_primers(forward_primers)
    filtered_reverse_primers = filter_primers(reverse_primers)

    primer_pairs = list(itertools.product(filtered_forward_primers, filtered_reverse_primers))

    # Filter out all pairs that have Tm's that are too far apart
    filtered_primer_pairs = [pair for pair in primer_pairs if math.isclose(pair[0].Tm, pair[1].Tm, abs_tol=args.D)]

    # Creates elements directly stored in output dict
    flat_primer_pair_info = []
    for i in filtered_primer_pairs:
        flat_primer_pair_info.append(
            {
                "Forward primer": i[0].get_primer_elements(),
                "Reverse primer": i[1].get_primer_elements()
            }
        )

    output_dict = {
        "Target gene": name_tag,
        "Genetic sequence": seq,
        "Total primer pairs": len(flat_primer_pair_info),
        "Generated with": "$software name",
        "Date": str(datetime.datetime.now()),
        "Primer Pair Info": flat_primer_pair_info
    }

    with open(args.o, 'w') as f:
        json.dump(output_dict, f, indent=2)

    print(f'The output file has been generated at {args.o}')


if __name__ == '__main__':
    main()
