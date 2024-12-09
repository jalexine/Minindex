import os
import sys
import pickle
import argparse
from collections import defaultdict
from time import time
from fmi import FmIndex

class Timer:
    """
    Timer class to measure elapsed time in seconds.
    """
    def __enter__(self):
        self.start = time()
        return self

    def __exit__(self, *args):
        self.end = time()
        self.elapsed = self.end - self.start

    def print(self, message):
        if hasattr(self, 'elapsed'):
            print(message.format(self.elapsed))
        else:
            print("Timer error: 'elapsed' attribute not set.")

def reverse_complement(kmer):
    """
    Returns the reverse complement of a k-mer.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement[base] for base in reversed(kmer))

def canonical_kmer(kmer):
    """
    Returns the canonical form of a k-mer (the lexicographically smallest between the k-mer and its reverse complement).
    """
    rev_comp = reverse_complement(kmer)
    return min(kmer, rev_comp)

def count_kmers(fasta_file, k):
    """
    Counts canonical k-mers in a FASTA file.
    """
    kmer_counts = defaultdict(int)
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                continue  # Skip headers
            sequence = line.strip().upper()
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                canonical = canonical_kmer(kmer)
                kmer_counts[canonical] += 1
    return kmer_counts

def filter_kmers(kmer_counts, threshold):
    """
    Filters k-mers with abundance below a given threshold.
    """
    return {kmer: count for kmer, count in kmer_counts.items() if count >= threshold}

def build_dbg(kmer_counts):
    """
    Builds an implicit de Bruijn graph from the k-mers.
    """
    dbg = defaultdict(list)
    for kmer in kmer_counts:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        dbg[prefix].append(suffix)
    return dbg

def generate_simplitigs(dbg):
    """
    Generates an SPSS using simplitigs from the de Bruijn graph.
    """
    visited = set()
    simplitigs = []

    for node in dbg:
        if node not in visited:
            simplitig = []
            stack = [node]

            # Extend forward
            while stack:
                current = stack.pop()
                if current not in visited:
                    simplitig.append(current[-1])
                    visited.add(current)
                    if current in dbg and len(dbg[current]) == 1:
                        stack.append(dbg[current][0])

            simplitigs.append("".join(simplitig))

    return "%".join(simplitigs) + "$"

def save_spss(spss_str, output_path, k, threshold, dataset_name):
    """
    Saves the SPSS string to a file.
    """
    filename = f"spss_{dataset_name}_k{k}_t{threshold}.pkl"
    filepath = os.path.join(output_path, filename)
    with open(filepath, "wb") as f:
        pickle.dump(spss_str, f)
    print(f"SPSS saved to {filepath}")

def save_fm_index(fm_index, output_path, k, threshold, dataset_name):
    """
    Saves the FM-index to a file.
    """
    filename = f"fmi_{dataset_name}_k{k}_t{threshold}.dump"
    filepath = os.path.join(output_path, filename)
    with open(filepath, "wb") as f:
        pickle.dump(fm_index, f)
    print(f"FM-index saved to {filepath}")

def save_statistics_and_print(output_path, fasta_file, k, t, time_selecting, time_spss, spss_size, num_spss, time_fmi):
    """
    Saves the statistics to a text file and prints them.
    """
    filename = os.path.join(output_path, f"stats_{os.path.basename(fasta_file)}.txt")
    content = (
        f"----------------------------------------\n"
        f"Processing: -i {fasta_file} -k {k} -t {t}\n"
        f"OUT TIME_SELECTING_KMERS={time_selecting:.2f} seconds\n"
        f"OUT TIME_SPSS_CONSTRUCTION={time_spss:.2f} seconds\n"
        f"OUT |SPSS(K)|={spss_size}\n"
        f"OUT #SPSS(K)={num_spss}\n"
        f"OUT TIME BUILD FMI={time_fmi:.2f} seconds\n"
        f"----------------------------------------\n"
    )
    
    # Print to terminal
    print(content)
    
    # Save to file
    with open(filename, "a") as f:
        f.write(content)

def test_fm_index(fm_index, spss):
    """
    Tests the FM-index by comparing with the raw SPSS string.
    """
    test_kmers = ["ACTG", "GGTA", "CCGT"]  # Example kmers for testing
    for kmer in test_kmers:
        in_fmi = fm_index.contains(kmer)
        in_raw = kmer in spss
        assert in_fmi == in_raw, f"FM-index inconsistency for {kmer}"
    print("FM-index validation passed.")

def main():
    parser = argparse.ArgumentParser(
        description="Generate SPSS from sequencing data and build FM-index.",
        add_help=False
    )
    parser.add_argument("-i", "--sequence", required=True, help="Sequence file (FASTA)")
    parser.add_argument("-k", "--kmer", type=int, required=True, help="K-mer size")
    parser.add_argument("-t", "--threshold", type=int, required=True, help="Threshold size") 
    parser.add_argument("-o", "--output", required=False, default="./output", 
                        help="Output directory for results (default: './output')")

    if len(sys.argv) == 1:
        print("\033[95m♡ pls use : sequences_to_indexed_spss.py -i SEQUENCE -k KMER -o OUTPUT [-t THRESHOLD]\033[0m")
        sys.exit(1)

    args = parser.parse_args()

    if not os.path.isfile(args.sequence):
        print(f"Error: Input sequence file '{args.sequence}' not found.")
        sys.exit(1)

    if not os.path.exists(args.output):
        print(f"Warning: Output directory '{args.output}' does not exist. Creating it.")
        os.makedirs(args.output)

    with Timer() as timer_select:
        kmer_counts = count_kmers(args.sequence, args.kmer)
        kmer_counts = filter_kmers(kmer_counts, args.threshold)
    
    dbg = build_dbg(kmer_counts)

    with Timer() as timer_spss:
        spss = generate_simplitigs(dbg)

    dataset_name = os.path.basename(args.sequence).split('.')[0]
    save_spss(spss, args.output, args.kmer, args.threshold, dataset_name)

    with Timer() as timer_build_fmi:
        fm_index = FmIndex(spss)
    save_fm_index(fm_index, args.output, args.kmer, args.threshold, dataset_name)

    # Validate FM-index
    test_fm_index(fm_index, spss)

    save_statistics_and_print(
        args.output,
        args.sequence,
        args.kmer,
        args.threshold,
        timer_select.elapsed,
        timer_spss.elapsed,
        len(spss),
        spss.count('%'),
        timer_build_fmi.elapsed,
    )

if __name__ == "__main__":
    main()
