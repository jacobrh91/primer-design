#!/usr/bin/env python3

#  Written by: Hans Müller Paul

# IMPORT LIBRARIES

import datetime
import itertools
import argparse
import math
import statistics
from os import path
from tabulate import tabulate ### only required for the table in output, can be replaced by something else

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
parser.add_argument('-m', '--mintemp', metavar='', dest='m', type=float, default=55,
                    help='min Tm in celsius, default = 55'
                    )
parser.add_argument('-x', '--maxtemp', metavar='', dest='x', type=float, default=62,
                    help='max Tm in celsius, default = 62'
                    )
parser.add_argument('-M', '--mingc', metavar='', dest='M', type=float, default=40,
                    help='min GC percentage, default = 40'
                    )
parser.add_argument('-X', '--maxgc', metavar='', dest='X', type=float, default=60,
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

def main():

    if args.verbose:
        print(args)

    # IMPORTING FASTA FILE

    if args.verbose:
        print('Item exists: ' + str(path.exists(args.i)))

    if path.exists(args.i):
        with open(args.i, 'r') as f:
            contents = f.read()
            seqlist = contents.split()
        if args.verbose:
            print('seqlist is {}'.format(type(seqlist)))        # Not sure this is valuable output, even in verbose mode
        name_tag = seqlist[0]
        seq = ''.join(seqlist[1:])
        if args.verbose:
            print(f'file name is {name_tag} and sequence is {seq}')
    else:
        print('file does not exist')

    # DEFINING REVERSE COMPLIMENT
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def revcomp(input):
        for k,v in alt_map.items():
            compseq = seq.replace(k,v)
        bases = list(compseq)
        bases = reversed([complement.get(base, base) for base in bases])
        bases = ''.join(bases)
        for k, v in alt_map.items():
            bases = bases.replace(v, k)
        return bases

    # DEFINING PRIMER LISTS

    def get_primer(input_string):
        length = len(input_string[:args.e])
        prim_list_1 = [input_string[i:j + 1] for i in range(length) for j in range(i+args.s,i+args.l)]
        if args.verbose:
            print('prim_list_1 is a {}'.format(type(prim_list_1)))
        return prim_list_1

    if args.verbose:
        print(get_primer(seq))

    def fwd(input):
        return get_primer(input)

    def rvs(input):
        seq_rc = revcomp(input)
        return get_primer(seq_rc)

    fwd_primers = fwd(seq)
    rvs_primers = rvs(seq)

    if args.verbose:
        print(f'List of forward primers: {fwd_primers}')
        print(f'List of reverse primers: {rvs_primers}')

    # DEFINING Primer GC CONTENT
    def GC(input):
        gc_count = input.count('g') + input.count('c') + input.count('G') + input.count('C')
        gc_fraction = float(gc_count)/len(input)
        return 100 * gc_fraction

    # DEFINING PRIMER MELTING TEMPERATURES
    def TM(input):
        N = len(input)
        melt = 81.5 + 0.41*(GC(input)) - (675/N)
        return melt

### NOTE FOR JACOB:
### The following two blocks were done in a very simplistic manner by removing sequences from a list
### if they did not meet a certain condition. I have done it in 4 separate steps: two for GC and two for
### TM, the first being less than minimum and the second being more than maximum. When I tried to implement
### a single step for verification (within the range) it returned an empty list. I may be attempting to do
### it in the wrong way.


    # DELETING FORWARD PRIMERS THAT DO NOT FALL WITHIN %GC AND TM LIMITS
    if args.verbose:
        print(f'Number of forward primers prior to GC verification: {len(fwd_primers)}')
    for sequence in fwd_primers:
        if GC(sequence) < args.M:
            fwd_primers.remove(sequence)
        elif GC(sequence) > args.X:
            fwd_primers.remove(sequence)

    if args.verbose:
        print(f'Number of forward primers after to GC verification: {len(fwd_primers)}')
        print(f'Number of forward primers prior to TM verification: {len(fwd_primers)}')
    for sequence in fwd_primers:
        if TM(sequence) < args.m:
            fwd_primers.remove(sequence)
        elif TM(sequence) > args.x:
            fwd_primers.remove(sequence)
    if args.verbose:
        print(f'Number of forward primers after to TM verification: {len(fwd_primers)}')

    # DELETING REVERSE PRIMERS THAT DO NOT FALL WITHIN %GC AND TM LIMITS

    if args.verbose:
        print(f'Number of reverse primers prior to GC verification: {len(rvs_primers)}')

    for sequence in rvs_primers:
        if GC(sequence) < args.M:
            rvs_primers.remove(sequence)
        elif GC(sequence) > args.X:
            rvs_primers.remove(sequence)
    if args.verbose:
        print(f'Number of reverse primers after to GC verification: {len(rvs_primers)}')
        print(f'Number of reverse primers prior to TM verification: {len(rvs_primers)}')

    for sequence in rvs_primers:
        if TM(sequence) < args.m    :
            rvs_primers.remove(sequence)
        elif TM(sequence) > args    .x:
            rvs_primers.remove(sequence)
    if args.verbose:
        print(f'Number of reverse primers after to TM verification: {len(rvs_primers)}')

    # FORMING PRIMER PAIRS

    if args.verbose:
        print(f"Pairing primers based on TM (primers in a pair will have a maximum TM difference of {args.D} degrees)")

    PrimerPairs = list(itertools.product(fwd_primers,rvs_primers))

    if args.verbose:
        print(f"Total number of primer pairs: {len(PrimerPairs)}")

    Pairs = [i for i in PrimerPairs if math.isclose(TM(i[0]),TM(i[1]),abs_tol=args.D)]

    if args.verbose:

        print(f"Number of primer pairs with matching TM: {len(Pairs)}")

    Temps = [(TM(i[0]),TM(i[1])) for i in Pairs]

    if args.verbose:

        print(f"TMs for each primer pair: {Temps}")

    MatchingTemp = list(zip(Pairs, Temps))

    if args.verbose:

        print(f"matched pairs: {MatchingTemp}")

    primertable = tabulate(MatchingTemp, headers=['Primer pair','TM of primer pair'])

    # WRITE TO OUTPUT FILE

    with open(args.o, 'w') as f2:
        f2.write(
            f"""
            Target gene: {name_tag}
            Genetic sequence:
            {seq}
            List of primers with TMs:
            {primertable}
            Generated with $software_name in {datetime.datetime.now()}
            Please cite $citation_quote
            """
        )

# CONFIRMATION MESSAGE

    print(f'The output file has been generated at {args.o}')


if __name__ == '__main__':
    main()
