import os
import pickle
import argparse
from collections import Counter
from timer import Timer
from fmi import FmIndex

# Pre-computed dictionary of complements
COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def reverse_complement(kmer):
    """Computes the reverse complement of a given k-mer."""
    return ''.join(COMPLEMENTS[nuc] for nuc in reversed(kmer))

def canonical_kmer(kmer):
    """Computes the canonical form of a k-mer."""
    rkmer = reverse_complement(kmer)
    return kmer if kmer < rkmer else rkmer

def process_sequence(sequence, counts, k):
    """Processes a sequence and updates the k-mers in the counter."""
    for i in range(len(sequence) - k + 1):
        kmer = canonical_kmer(sequence[i:i + k])
        counts[kmer] += 1

def kmer_abundance(fasta_file, k):
    """
    Computes the abundance of all canonical k-mers of size k from a FASTA file.
    """
    counts = Counter()
    with open(fasta_file, "r") as f:
        sequence = ""
        for line in f:
            if line.startswith(">"):  
                if sequence:
                    process_sequence(sequence, counts, k)
                sequence = ""
            else:
                sequence += line.strip().upper()
        if sequence:
            process_sequence(sequence, counts, k)
    return counts

def neighbors(kmer, dbg_set, is_right):
    """
    Generates neighbors to the left or right of a k-mer.
    Explores all possible paths to maximize SPSS fragments.
    """
    base_kmer = kmer[1:] if is_right else kmer[:-1]
    return [base_kmer + c if is_right else c + base_kmer
            for c in "ACGT" if (base_kmer + c if is_right else c + base_kmer) in dbg_set]

def create_spss(dbg):
    """
    Builds SPSS by maximizing the diversity of explored paths.
    """
    spss_fragments = []
    dbg_set = set(dbg.keys())

    while dbg_set:
        kmer = dbg_set.pop()
        left_part, right_part = [], []
        current = kmer

        # Right extension
        while True:
            neigh = neighbors(current, dbg_set, True)
            if neigh:
                # Explore all neighbors to maximize diversity
                next_kmer = max(neigh, key=lambda n: dbg.get(n, 0))  # Choose the most frequent neighbor
                right_part.append(next_kmer[-1])
                dbg_set.remove(next_kmer)
                current = next_kmer
            else:
                break

        current = kmer
        # Left extension
        while True:
            neigh = neighbors(current, dbg_set, False)
            if neigh:
                # Explore all neighbors to maximize diversity
                prev_kmer = max(neigh, key=lambda n: dbg.get(n, 0))  # Choose the most frequent neighbor
                left_part.append(prev_kmer[0])
                dbg_set.remove(prev_kmer)
                current = prev_kmer
            else:
                break

        # Build the unitig
        unitig = ''.join(reversed(left_part)) + kmer + ''.join(right_part)
        spss_fragments.append(unitig)

    return '%'.join(spss_fragments) + '$'

def save_spss(spss_str, output_path):
    """
    Saves the SPSS string to a file using pickle.
    """
    filepath = os.path.join(output_path, "SPSS_seq.pkl") if os.path.isdir(output_path) else output_path
    with open(filepath, "wb") as f:
        pickle.dump(spss_str, f)

def main():
    parser = argparse.ArgumentParser(description="Generate SPSS from sequencing data.")
    parser.add_argument("-i", "--sequence", required=True, help="Sequence file (FASTA)")
    parser.add_argument("-t", "--threshold", type=int, required=False, default=1,
                        help="Minimum abundance threshold (not used in this optimized version)")
    parser.add_argument("-k", "--kmer", type=int, required=True, help="K-mer size")
    parser.add_argument("-o", "--output", required=True, help="Output file or directory")
    args = parser.parse_args()

    # Count k-mers
    with Timer() as timer_select:
        dbg = kmer_abundance(args.sequence, args.kmer)
        print(f"Total number of k-mers: {len(dbg)}")
    timer_select.print("OUT TIME_SELECTING_KMERS={} seconds")

    # Build SPSS
    with Timer() as timer_spss:
        spss = create_spss(dbg)
        print(f"Size of SPSS string: {len(spss)} characters")
    timer_spss.print("OUT TIME_SPSS_CONSTRUCTION={} seconds")

    # Save SPSS
    save_spss(spss, args.output)

    # Statistics
    total_characters = len(spss)
    distinct_sequences = spss.count('%')
    print(f"OUT |SPSS(K)|={total_characters}")
    print(f"OUT #SPSS(K)={distinct_sequences}")

    # Build FM-Index
    with Timer() as timer_build_fmi:
        fm_index = FmIndex(spss)
        fm_index.save(args.output)
    timer_build_fmi.print("OUT TIME BUILD FMI={} seconds")

if __name__ == "__main__":
    main()
