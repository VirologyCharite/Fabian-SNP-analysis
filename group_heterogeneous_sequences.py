#!/usr/bin/env python3

from dark.reads import addFASTACommandLineOptions, parseFASTACommandLineOptions
from dark.fasta import FastaReads
from math import log10
from heterogeneous_sequences import heterogeneousSites
import argparse
import sys

parser = argparse.ArgumentParser()

parser.add_argument("--homogenFraction", type=float, default=1.0,
                    help="If the fraction of the most common nucleotide at a"
                    "site is bigger than this value, the site will be considered homogenous.")

parser.add_argument("--bestSNP", default=False, action='store_true',
                    help="Returns the SNP where the homogenFraction is smallest.")

addFASTACommandLineOptions(parser)

args = parser.parse_args()

homogenFraction = args.homogenFraction
bestSNP = args.bestSNP
reads = list(parseFASTACommandLineOptions(args))

if not 0 < homogenFraction <= 1:
    raise ValueError('--homogenFraction needs to be between 0 and 1.')

if len(set(len(read) for read in reads)) != 1:
    raise ValueError('Not all read lengths are the same.....')
else:
    length = len(reads[0])

    count, ids, indexes = heterogeneousSites(reads, length, homogenFraction)

    width = int(log10(length)) + 1
    for index in indexes:
        print('A SNP is detected at index %*d.' % (width, index + 1))
        print('Detected bases: %s' % (str(count[index])[9:-2]))

        if 'A' in count[index].keys():
            print('A:')
            for id in ids[index]['A']:
                print(id)

        if 'C' in count[index].keys():
            print('C:')
            for id in ids[index]['C']:
                print(id)

        if 'G' in count[index].keys():
            print('G:')
            for id in ids[index]['G']:
                print(id)

        if 'T' in count[index].keys():
            print('T:')
            for id in ids[index]['T']:
                print(id)
if bestSNP:
    print('The SNP that where the nucleotides are 50:50.')

