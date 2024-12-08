#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FM-index implementation for genome indexing and querying.
"""

import pickle
from pysuffix3 import tools_karkkainen_sanders as tks

def load_fm_index(filename):
    """
    Load the FM-index from a file and return it.
    """
    with open(filename, 'rb') as f:
        fm_index = pickle.load(f)
    return fm_index

class FmIndex:
    """
    Class for storing an FM-index of a sequence.
    """
    def __init__(self, sequence):
        """
        Constructor: Initializes the FM-index.
        """
        # Directly use tks to generate the suffix array
        self.suffix_array = tks.simple_kark_sort(sequence)
        #print("Suffix Array:", self.suffix_array)

        self.sequence = sequence
        self.bwt = self.set_bwt()
        self.rank, self.n = self.set_n_and_rank()
        #print("N:", self.n)
        #print("Rank:", self.rank)

    def save(self, filename):
        """
        Save the FM-index to a file using pickle.
        """
        with open(filename, 'wb') as f:
            pickle.dump(self, f)
        #print(f"FM-index saved to {filename}.")

    def set_bwt(self):
        """
        Given a sequence and its suffix array, computes the Burrows-Wheeler Transform (BWT).
        """
        bwt = "".join(self.sequence[i - 1] if i > 0 else self.sequence[-1] for i in self.suffix_array)
        #print("Corrected BWT:", bwt)
        return bwt

    def set_n_and_rank(self):
        """
        Computes N (occurrences of each letter in the BWT) and rank (cumulative counts for each character).
        """
        bw = self.bwt
        tots = dict()
        ranks = []

        for c in bw:
            if c not in tots:
                tots[c] = 0
            ranks.append(tots[c])
            tots[c] += 1

        #print("N (Total counts):", tots)
        #print("Ranks:", ranks)
        return ranks, tots

    def lf(self, alpha, k) -> int:
        """
        Returns the index in the suffix array corresponding to the k-th suffix starting with letter alpha.
        """
        lf_index = self.n[alpha] + k - 1
        #print(f"LF Mapping for ({alpha}, {k}): {lf_index}")
        return lf_index

    def find_next(self, alpha, l) -> int:
        """
        Find the first line >= l such that BWT[line] == alpha.
        """
        for i in range(l, len(self.bwt)):
            if self.bwt[i] == alpha:
                #print(f"Find Next: Found {alpha} at index {i}")
                return i
        #print(f"Find Next: {alpha} not found after index {l}")
        return -1

    def find_prev(self, alpha, l) -> int:
        """
        Find the last line <= l such that BWT[line] == alpha.
        """
        for i in range(l, -1, -1):
            if self.bwt[i] == alpha:
                #print(f"Find Prev: Found {alpha} at index {i}")
                return i
        #print(f"Find Prev: {alpha} not found before or at index {l}")
        return -1
    
    def contains(self, q) -> bool:
        """
        Check if the query q is indexed in the FM-index.
        """
        # Sort characters and build the C array
        chars = sorted(self.n.keys())
        C = {}
        total = 0
        for c in chars:
            C[c] = total
            total += self.n[c]

        # Helper to compute occurrences of `c` up to index `i`
        def occ(c, i):
            if i < 0:
                return 0
            return sum(1 for j in range(i + 1) if self.bwt[j] == c)

        # Initialize the range [l, r]
        l, r = 0, len(self.bwt) - 1

        # Process the query in reverse
        for char in reversed(q):
            if char not in C:
                #print(f"Character {char} not in BWT. Query not found.")
                return False
            l = C[char] + occ(char, l - 1)
            r = C[char] + occ(char, r) - 1
            #print(f"Updated range for {char}: l={l}, r={r}")
            if l > r:
                #print(f"Range invalid after processing {char}. Query not found.")
                return False

        return l <= r

'''
def main():
    sequence = "banana$"
    fm_index = FmIndex(sequence)

    fm_index.save("test_fm_index.dump")

    print("Final BWT:", fm_index.bwt)
    print("Final N:", fm_index.n)
    print("Final Ranks:", fm_index.rank)

    # Test contains
    print("Contains Tests:")
    print("Contains 'ana':", fm_index.contains('ana'))
    print("Contains 'banana':", fm_index.contains('banana'))
    print("Contains 'nana':", fm_index.contains('nana'))
    print("Contains 'xyz':", fm_index.contains('xyz'))


if __name__ == "__main__":
    main()
'''
