#!/usr/bin/env python3

from __future__ import print_function, division
from collections import Counter, defaultdict
from dark.fasta import FastaReads

def heterogeneousSites(reads, length, homogeneous_fraction):
    """
    Returns the sites of heterogeneity in a list of FastaReads of the same length.

    @param reads: An iterable of C{dark.reads.Read} instances.
    @param length: The length of each read in reads.
    @param homogeneous_fraction: Cutoff frequency at which a site is considered to be homogeneous.
    @raises NameError: If C{reads} is empty.
    @return: A C{dict} with C{int} index keys and C{Counter} instances
        as values.
    """
    result = {}
    resultids = {}
    heterogeneousIndexes = []

    for index in range(length):
        counts = Counter()
        ids = {}

        for read in reads:
            # Ignore gaps.
            if read.sequence[index] != '-':
                counts[read.sequence[index]] += 1
                # give ids as output
                if read.sequence[index] in ids:
                    ids[read.sequence[index]].append(read.id)
                if read.sequence[index] not in ids:
                    ids[read.sequence[index]] = [read.id]

        if len(counts) > 1 and counts.most_common(1)[0][1] / sum(counts.values()) <= homogeneous_fraction:
            result[index] = counts
            resultids[index] = ids
            heterogeneousIndexes.append(index)

    return result, resultids, heterogeneousIndexes
