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
    print("LF('$', 1):", fm_index.lf('$', 1))  # Résultat attendu: 0 (n['$']=1, 1+1-1=1 mais comme c'est 0-based, à vérifier)
    # Note : Ici, il faudra vérifier ce résultat. Selon la logique du code, lf('$',1)=1+1-1=1. Si on s'attend à un index 0-based
    # il peut y avoir confusion. On laisse le résultat brut du code tel quel.

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

    # Test de la sérialisation et désérialisation
    fm_index.save("banana_fm_index.dump")
    loaded_fm = load_fm_index("banana_fm_index.dump")
    print("Loaded FM-index BWT:", loaded_fm.bwt)
    # Résultat attendu: la même BWT que fm_index.bwt, soit "annb$aa"



if __name__ == "__main__":
    main()

