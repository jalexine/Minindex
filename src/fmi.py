#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
FM-index implementation for genome indexing and querying.
"""

import pickle
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), "../pysuffix3"))
from pysuffix3 import tools_karkkainen_sanders as tks

def load_fm_index(filename):
    """
    Load an FM-index object from a given file.

    Parameters:
        filename (str): The path to the file containing the serialized FM-index object.

    Returns:
        FmIndex: The deserialized FM-index instance loaded from the file.
    """
    with open(filename, 'rb') as f:
        fm_index = pickle.load(f)
    return fm_index

class FmIndex:
    """
    A class that constructs and stores an FM-index for efficient substring queries.

    Parameters:
        sequence (str): The input sequence (e.g., a genome) for which the FM-index will be built.

    Attributes:
        suffix_array (list[int]): The suffix array of the input sequence.
        sequence (str): The original input sequence.
        bwt (str): The Burrows–Wheeler transform of the input sequence.
        rank (list[int]): The array of rank values associated with each character in the BWT.
        n (dict): A dictionary mapping each character to its total count in the BWT.
    """

    def __init__(self, sequence):
        """
        Initialize the FM-index by computing the suffix array, BWT, and related structures.

        Parameters:
            sequence (str): The input sequence to index.
        """
        self.suffix_array = tks.simple_kark_sort(sequence)
        self.sequence = sequence
        self.bwt = self.set_bwt()
        self.rank, self.n = self.set_n_and_rank()

    def save(self, filename):
        """
        Serialize the current FM-index instance and save it to a file.

        Parameters:
            filename (str): The path of the file where the FM-index should be saved.

        Returns:
            None
        """
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    def set_bwt(self):
        """
        Compute the Burrows–Wheeler transform (BWT) from the suffix array and sequence.

        Parameters:
            None

        Returns:
            str: The BWT string derived from the input sequence.
        """
        bwt = "".join(self.sequence[i - 1] if i > 0 else self.sequence[-1] for i in self.suffix_array)
        return bwt

    def set_n_and_rank(self):
        """
        Compute the counts of each character and rank arrays from the BWT.

        Parameters:
            None

        Returns:
            tuple:
                - list[int]: The rank array of characters in the BWT.
                - dict: A dictionary where each key is a character and its value is the count of that character in the BWT.
        """
        bw = self.bwt
        tots = {}
        ranks = []

        for c in bw:
            if c not in tots:
                tots[c] = 0
            ranks.append(tots[c])
            tots[c] += 1

        return ranks, tots

    def lf(self, alpha, k):
        """
        Compute the LF-mapping for a given character and occurrence count (k).

        LF-mapping finds the position in the suffix array that corresponds to the k-th occurrence of character alpha in the BWT.

        Parameters:
            alpha (str): The character to map.
            k (int): The occurrence count of the character alpha.

        Returns:
            int: The index in the suffix array corresponding to the given character occurrence.
        """
        lf_index = self.n[alpha] + k - 1
        return lf_index

    def find_next(self, alpha, l):
        """
        Find the next occurrence of a given character in the BWT starting from index l.

        Parameters:
            alpha (str): The character to find in the BWT.
            l (int): The starting index in the BWT from which to begin the search.

        Returns:
            int: The index of the next occurrence of alpha in the BWT at or after l, or -1 if not found.
        """
        for i in range(l, len(self.bwt)):
            if self.bwt[i] == alpha:
                return i
        return -1

    def find_prev(self, alpha, l):
        """
        Find the previous occurrence of a given character in the BWT up to index l.

        Parameters:
            alpha (str): The character to find in the BWT.
            l (int): The ending index in the BWT at which to end the backward search.

        Returns:
            int: The index of the last occurrence of alpha at or before l, or -1 if not found.
        """
        for i in range(l, -1, -1):
            if self.bwt[i] == alpha:
                return i
        return -1
    
    def occ(self, c, i):
        """
        Count the number of occurrences of character c in the BWT from index 0 to index i inclusive.

        Parameters:
            c (str): The character whose occurrences are to be counted.
            i (int): The ending index of the range in the BWT (0-based).

        Returns:
            int: The number of occurrences of c in the BWT up to and including index i.
        """
        if i < 0:
            return 0
        return sum(1 for j in range(i + 1) if self.bwt[j] == c)

    def contains(self, q):
        """
        Check whether a given query substring q exists in the indexed sequence.

        The method uses backward search on the BWT to determine if q occurs in the sequence.

        Parameters:
            q (str): The query substring to search for.

        Returns:
            bool: True if q is found in the sequence, False otherwise.
        """
        # Sort characters and build the C array (cumulative counts)
        chars = sorted(self.n.keys())
        C = {}
        total = 0
        for c in chars:
            C[c] = total
            total += self.n[c]

        # Initialize the search range [l, r]
        l, r = 0, len(self.bwt) - 1

        # Process the query in reverse
        for char in reversed(q):
            if char not in C:
                return False
            l = C[char] + self.occ(char, l - 1)
            r = C[char] + self.occ(char, r) - 1
            if l > r:
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

    # Tests supplémentaires
    print("Contains 'ban':", fm_index.contains('ban'))     # True, 'ban' est le début de "banana$"
    print("Contains 'na':", fm_index.contains('na'))       # True, 'na' apparaît plusieurs fois, ex: "ba(n a)na$"
    print("Contains 'aa':", fm_index.contains('aa'))       # False, il n'y a jamais deux 'a' consécutifs
    print("Contains '$':", fm_index.contains('$'))         # True, '$' est le symbole de fin dans "banana$"
    print("Contains 'na$':", fm_index.contains('na$'))     # True, 'na$' se trouve à la fin: "bana(n a $)"
    print("Contains 'banan':", fm_index.contains('banan')) # True, 'banan' est le début de la séquence "banana"
    print("Contains 'anana':", fm_index.contains('anana')) # True, 'anana' est présente: "b(anana)$"
    print("Contains 'bana$':", fm_index.contains('bana$')) # False, il n'y a pas de sous-chaîne 'bana$' en continu

        # Tests pour les autres méthodes de la classe

    # Test de la fonction lf (LF-mapping)
    # Avec la séquence "banana$", la BWT devrait ressembler à quelque chose comme "annb$aa" (selon le tri).
    # Comptage des caractères dans la BWT : a:3, n:2, b:1, $:1
    # lf(alpha, k) = n[alpha] + k - 1
    print("LF('a', 1):", fm_index.lf('a', 1))  # Résultat attendu: 3 (puisque n['a']=3, 3+1-1=3)
    print("LF('a', 2):", fm_index.lf('a', 2))  # Résultat attendu: 4 (3+2-1=4)
    print("LF('n', 1):", fm_index.lf('n', 1))  # Résultat attendu: 2 (n['n']=2, 2+1-1=2)
    print("LF('n', 2):", fm_index.lf('n', 2))  # Résultat attendu: 3 (2+2-1=3)
    print("LF('$', 1):", fm_index.lf('$', 1))  # Résultat attendu: 1 (n['$']=1, 1+1-1=1 mais comme c'est 0-based, à vérifier)

    # Test de find_next et find_prev
    # BWT: 'annb$aa' (index: a(0), n(1), n(2), b(3), $(4), a(5), a(6))
    print("find_next('a', 0):", fm_index.find_next('a', 0))  # Cherche 'a' à partir de l'indice 0, trouvé à 0
    # Résultat attendu: 0
    print("find_next('b', 0):", fm_index.find_next('b', 0))  # Cherche 'b' depuis 0, b est à l'indice 3
    # Résultat attendu: 3
    print("find_next('$', 2):", fm_index.find_next('$', 2))   # Cherche '$' à partir de 2, trouvé à 4
    # Résultat attendu: 4

    print("find_prev('a', 6):", fm_index.find_prev('a', 6))  # Cherche 'a' avant ou à 6, trouvé à 6
    # Résultat attendu: 6
    print("find_prev('b', 6):", fm_index.find_prev('b', 6))  # Cherche 'b' avant ou à 6 (donc 6,5,4,3...), trouvé à 3
    # Résultat attendu: 3
    print("find_prev('$', 2):", fm_index.find_prev('$', 2))   # Cherche '$' avant ou à 2, indices 2('n'),1('n'),0('a'), pas de '$'
    # Résultat attendu: -1 (non trouvé)




if __name__ == "__main__":
    main()
'''
