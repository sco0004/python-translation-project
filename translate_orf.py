#!/bin/env python3

import sys
import re

def main():
    import argparse
    import find_orf
    import translate

# Creating the command-line parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

# Setting up what command-line arguments the parser can receive
    parser.add_argument('-sequence',
            metavar= 'SEQUENCE',
            type=str,
            help=('The sequence to search for an open-reading frame. If the path flag ('-p'/'--path') is specified, then this should be a path to a file containing the sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = (('The sequence argument should be treated as a path to a containing the sequence to be searched. default: False'))
    parser.add_argument('-s', '--START_CODONS', '--start_codons',
            action = 'store_true',
            help = (('One or more possible start codons. (default: [AUG])'))

    parser.add_argument('-x', '--STOP_CODONS', '--stop_codons',
            action = 'store_true',
            help = (('One or more possible stop codons. (default: [UAA,UAG, UGA]'))

# Parse the comman-line arguments
    args = parser.parse_args()

    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    orf = find_orf.find_first_orf(sequence = sequence,
            start_codons = args.start_codons,
            stop_codons = args.stop_codons)

    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}

    translation = translate.translate_sequence(rna_sequence = orf, genetic_code = genetic_code)
    sys.stdout.write('{}\n'.format(translation))

if __name__ == '__main__':
    main()
