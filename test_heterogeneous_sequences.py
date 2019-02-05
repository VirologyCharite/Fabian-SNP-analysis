#!/usr/bin/env python3

from unittest import TestCase
from dark.reads import Read, Reads
from heterogeneous_sequences import heterogeneousSites

class TestHeterogeneousSites(TestCase):
    """
    Test the heterogeneousSites function.
    """

    def testOneRead(self):
        """
        heterogeneousSites must return an empty dictionary if only one read is given.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read])

        self.assertEqual(({}, {}, []), heterogeneousSites(reads, len(read), 1))

    def testHomogeneousReads(self):
        """
        heterogeneousSites must return an empty dictionary if homogenous reads are given.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCG')])

        self.assertEqual(({}, {}, []), heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsOneDifference(self):
        """
        heterogeneousSites must return a dictionary with one entry as expected if 
        reads given differ at one site.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC')])

        self.assertEqual(({3: {'G': 1, 'C': 1}},
                          {3: {'G': ['id'], 'C': ['id2']}}, [3]),
                           heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsTwoDifferences(self):
        """
        heterogeneousSites must return a dictionary with two entries as expected if 
        reads given differ at two sites.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'TCCC')])

        self.assertEqual(({0: {'A': 1, 'T': 1}, 3: {'G': 1, 'C': 1}},
                          {0: {'A': ['id'], 'T': ['id2']},
                           3: {'C': ['id2'], 'G': ['id']}}, [0, 3]),
                           heterogeneousSites(reads, len(read), 1))

    def testHeterogeneousReadsFractionHigh(self):
        """
        heterogeneousSites must return a dictionary with one entry as expected if 
        reads given differ and are less homogeneous than specified by the homogeneity 
        cutoff fraction.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC'), Read('id3', 'ACCC')])

        self.assertEqual(({3: {'C': 2, 'G': 1}},
                          {3: {'G': ['id'], 'C': ['id2', 'id3']}}, [3]),
                           heterogeneousSites(reads, len(read), 0.7))

    def testHeterogeneousReadsFractionLow(self):
        """
        heterogeneousSites must return an empty dictionary if reads given differ and 
        are more homogeneous than specified by the homogeneity cutoff fraction.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'ACCC'), Read('id3', 'ACCC')])

        self.assertEqual(({}, {}, []), heterogeneousSites(reads, len(read), 0.6))

    def testHeterogeneousReadsFractionLowWithOneDifference(self):
        """
        heterogeneousSites must return a dictionary with one entry if reads given differ 
        at two sites and at one site are more homogeneous than specified by the homogeneity 
        cutoff fraction; at the other site less homogeneous than specified by the 
        homogeneity cutoff fraction.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', 'TCCG'), Read('id3', 'TCCG'), Read('id4', 'ACCG')])

        self.assertEqual(({0: {'A': 2, 'T': 2}},
                          {0: {'A': ['id', 'id4'], 'T': ['id2', 'id3']}}, [0]),
                           heterogeneousSites(reads, len(read), 0.6))

    def testGaps(self):
        """
        heterogeneousSites must return an empty dictionary if reads given differ only by gaps;
        gaps do not count towards heterogeneity.
        """
        read = Read('id', 'ACCG')
        reads = Reads([read, Read('id2', '-CC-')])

        self.assertEqual(({}, {}, []), heterogeneousSites(reads, len(read), 1))

    def testEmptyReads(self):
        """
        heterogeneousSites must return an empty dictionary if no reads are given.
        """
        reads = Reads()

        self.assertEqual(({}, {}, []), heterogeneousSites(reads, 0, 1))
