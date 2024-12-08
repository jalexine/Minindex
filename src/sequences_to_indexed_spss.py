import os
import sys
import pickle
import argparse
from collections import Counter
from time import time
from fmi import FmIndex

COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

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
            neighbors = [
                current[1:] + c for c in "ACGT" if current[1:] + c in dbg_set
            ]
            if neighbors:
                next_kmer = max(neighbors, key=lambda n: dbg.get(n, 0))  # Most frequent neighbor
                right_part.append(next_kmer[-1])
                dbg_set.remove(next_kmer)
                current = next_kmer
            else:
                break

        current = kmer
        # Left extension
        while True:
            neighbors = [
                c + current[:-1] for c in "ACGT" if c + current[:-1] in dbg_set
            ]
            if neighbors:
                prev_kmer = max(neighbors, key=lambda n: dbg.get(n, 0))  # Most frequent neighbor
                left_part.append(prev_kmer[0])
                dbg_set.remove(prev_kmer)
                current = prev_kmer
            else:
                break

        # Build the unitig
        unitig = ''.join(reversed(left_part)) + kmer + ''.join(right_part)
        spss_fragments.append(unitig)

    return '%'.join(spss_fragments) + '$'

def save_spss(spss_str, output_path, k, threshold, dataset_name):
    """
    Saves the SPSS string to a file using pickle with an explicit name.
    """
    filename = f"spss_{dataset_name}_k{k}_t{threshold}.pkl"
    filepath = os.path.join(output_path, filename)
    with open(filepath, "wb") as f:
        pickle.dump(spss_str, f)
    print(f"SPSS saved to {filepath}")

def save_fm_index(fm_index, output_path, k, threshold, dataset_name):
    """
    Saves the FM-index to a file with a clear and explicit name.
    """
    filename = f"fmi_{dataset_name}_k{k}_t{threshold}.dump"
    filepath = os.path.join(output_path, filename)
    with open(filepath, "wb") as f:
        pickle.dump(fm_index, f)
    print(f"FM-index saved to {filename}")

def save_statistics_and_print(output_path, fasta_file, k, t, time_selecting, time_spss, spss_size, num_spss, time_fmi):
    """
    Saves the statistics to a text file and prints them to the terminal.
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

    # Custom help message
    if len(sys.argv) == 1:
        print("\033[95m♡ pls use : sequences_to_indexed_spss.py -i SEQUENCE -k KMER -o OUTPUT [-t THRESHOLD] <n>\033[0m")
        sys.exit(1)

    args = parser.parse_args()

    # Validate input file exists
    if not os.path.isfile(args.sequence):
        print(f"Error: Input sequence file '{args.sequence}' not found.")
        sys.exit(1)

    # Validate output directory
    if not os.path.exists(args.output):
        print(f"Warning: Output directory '{args.output}' does not exist. Creating it.")
        os.makedirs(args.output)

    # Core processing
    with Timer() as timer_select:
        dbg = kmer_abundance(args.sequence, args.kmer)
        dbg = {kmer: count for kmer, count in dbg.items() if count >= args.threshold}

    with Timer() as timer_spss:
        spss = create_spss(dbg)

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
