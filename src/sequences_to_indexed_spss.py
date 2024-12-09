import os
import sys
import pickle
import argparse
import csv
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
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join(complement[base] for base in reversed(kmer))

def canonical_kmer(kmer):
    rev_comp = reverse_complement(kmer)
    return min(kmer, rev_comp)

def count_kmers(fasta_file, k):
    kmer_counts = defaultdict(int)
    with open(fasta_file, "r") as file:
        for line in file:
            if line.startswith(">"):
                continue
            sequence = line.strip().upper()
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                canonical = canonical_kmer(kmer)
                kmer_counts[canonical] += 1
    return kmer_counts

def filter_kmers(kmer_counts, threshold):
    return {kmer: count for kmer, count in kmer_counts.items() if count >= threshold}

def build_dbg(kmer_counts):
    dbg = defaultdict(list)
    for kmer in kmer_counts:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        dbg[prefix].append(suffix)
    return dbg

def generate_simplitigs(dbg):
    visited = set()
    simplitigs = []

    for node in dbg:
        if node not in visited:
            simplitig = []
            stack = [node]

            while stack:
                current = stack.pop()
                if current not in visited:
                    simplitig.append(current[-1])
                    visited.add(current)
                    if current in dbg and len(dbg[current]) == 1:
                        stack.append(dbg[current][0])

            simplitigs.append("".join(simplitig))

    return "%".join(simplitigs) + "$"

def generate_unitigs(dbg):
    """
    Generates unitigs from the de Bruijn graph.
    """
    visited = set()
    unitigs = []

    nodes = list(dbg.keys())

    for node in nodes:
        if node not in visited:
            unitig = [node]
            visited.add(node)
            next_node = dbg[node][0] if len(dbg[node]) == 1 else None

            while next_node and next_node not in visited and len(dbg[next_node]) == 1:
                unitig.append(next_node[-1])
                visited.add(next_node)
                next_node = dbg[next_node][0]

            unitigs.append("".join(unitig))
    return "%".join(unitigs) + "$"


def save_fm_index(fm_index, output_base, k, threshold, dataset_name, mode):
    """
    Saves the FM-index to a file in the correct mode-specific directory.
    """
    output_path = output_base
    os.makedirs(output_path, exist_ok=True)  # Crée uniquement les dossiers nécessaires
    filename = f"fmi_{dataset_name}_k{k}_t{threshold}_{mode}.dump"
    filepath = os.path.join(output_path, filename)
    with open(filepath, "wb") as f:
        pickle.dump(fm_index, f)
    print(f"FM-index saved to {filepath}")


def save_statistics_and_print(output_base, fasta_file, k, t, time_selecting, time_spss, spss_size, num_spss, time_fmi, mode):
    """
    Saves the statistics to a CSV file in the correct mode-specific directory.
    """
    parent_directory = os.path.abspath(os.path.join(output_base, os.pardir))
    filename = os.path.join(parent_directory, f"stats_{os.path.basename(fasta_file)}.csv")
    headers = ["Dataset", "k", "t", "TIME_SELECTING_KMERS", "TIME_SPSS_CONSTRUCTION", "SPSS(K)", "#SPSS(K)", "TIME_BUILD_FMI"]
    data = [
        os.path.basename(fasta_file),
        k,
        t,
        f"{time_selecting:.2f}",
        f"{time_spss:.2f}",
        spss_size,
        num_spss,
        f"{time_fmi:.2f}",
    ]
    write_headers = not os.path.isfile(filename)
    
    with open(filename, mode='a', newline='') as file:
        writer = csv.writer(file)
        if write_headers:
            writer.writerow(headers)
        writer.writerow(data)

    content = (
        f"----------------------------------------\n"
        f"Processing: -i {fasta_file} -k {k} -t {t} --mode {mode}\n"
        f"OUT TIME_SELECTING_KMERS={time_selecting:.2f} seconds\n"
        f"OUT TIME_SPSS_CONSTRUCTION={time_spss:.2f} seconds\n"
        f"OUT |SPSS(K)|={spss_size}\n"
        f"OUT #SPSS(K)={num_spss}\n"
        f"OUT TIME BUILD FMI={time_fmi:.2f} seconds\n"
        f"----------------------------------------\n"
    )
    print(content)



def test_fm_index(fm_index, spss):
    test_kmers = ["ACTG", "GGTA", "CCGT"]
    for kmer in test_kmers:
        in_fmi = fm_index.contains(kmer)
        in_raw = kmer in spss
        assert in_fmi == in_raw, f"FM-index inconsistency for {kmer}"
    print("FM-index validation passed.")

def main():
    parser = argparse.ArgumentParser(description="Generate SPSS from sequencing data and build FM-index.")
    parser.add_argument("-i", "--sequence", required=True, help="Sequence file (FASTA)")
    parser.add_argument("-k", "--kmer", type=int, required=True, help="K-mer size")
    parser.add_argument("-t", "--threshold", type=int, required=True, help="Threshold size")
    parser.add_argument("-o", "--output", default="./Benchmark", help="Output directory (default: './Benchmark')")
    parser.add_argument("--mode", default="simplitigs", choices=["simplitigs", "unitigs"],
                        help="Mode of construction: 'simplitigs' (default) or 'unitigs'")

    args = parser.parse_args()

    if not os.path.isfile(args.sequence):
        print(f"Error: Input sequence file '{args.sequence}' not found.")
        sys.exit(1)

    with Timer() as timer_select:
        kmer_counts = count_kmers(args.sequence, args.kmer)
        kmer_counts = filter_kmers(kmer_counts, args.threshold)
    
    dbg = build_dbg(kmer_counts)

    with Timer() as timer_spss:
        if args.mode == "unitigs":
            spss = generate_unitigs(dbg)
        else:
            spss = generate_simplitigs(dbg)

    dataset_name = os.path.basename(args.sequence).split('.')[0]

    with Timer() as timer_build_fmi:
        fm_index = FmIndex(spss)
    save_fm_index(fm_index, args.output, args.kmer, args.threshold, dataset_name, args.mode)

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
        args.mode
    )


if __name__ == "__main__":
    main()