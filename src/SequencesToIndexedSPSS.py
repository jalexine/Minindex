from collections import defaultdict
from fmi import FmIndex

def reverse_complement(kmer):
    """
    Returns the reverse complement of a k-mer.
    
    Args:
        kmer (str): The input k-mer.
    
    Returns:
        str: The reverse complement of the k-mer.
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp = ''.join(complement[base] for base in reversed(kmer))
    return rev_comp

def canonical_kmer(kmer):
    """
    Returns the canonical form of a k-mer.
    The canonical form is the lexicographically smaller string 
    between the k-mer and its reverse complement.
    
    Args:
        kmer (str): The input k-mer.
    
    Returns:
        str: The canonical form of the k-mer.
    """
    rev_comp = reverse_complement(kmer)
    return min(kmer, rev_comp)

def count_solid_kmers(sequence, k, solidity_threshold):
    """
    Counts k-mers in a sequence and filters them based on a solidity threshold.
    Only canonical k-mers with occurrences >= solidity_threshold are retained.

    Args:
        sequence (str): The input DNA sequence.
        k (int): Length of the k-mers.
        solidity_threshold (int): Minimum number of occurrences to consider a k-mer solid.
    
    Returns:
        set: A set of solid canonical k-mers.
    """
    kmer_count = defaultdict(int)
    n = len(sequence)

    # Sliding window to count k-mers
    for i in range(k, n):
        kmer = kmer[1:] + sequence[i] 
        canonical = canonical_kmer(kmer)
        kmer_count[canonical] += 1

    # Filter k-mers by solidity threshold
    solid_kmers = {kmer for kmer, count in kmer_count.items() if count >= solidity_threshold}

    return solid_kmers

def build_de_bruijn_graph(solid_kmers):
    """
    Constructs a De Bruijn graph from a set of solid k-mers.
    
    Args:
        solid_kmers (set): A set of solid canonical k-mers.
    
    Returns:
        defaultdict: A graph represented as adjacency lists.
    """
    graph = defaultdict(list)
    for kmer in solid_kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
    return graph

def make_node_edge_map(graph):
    """
    Creates a mapping of graph nodes to their neighbors (edges).
    
    Args:
        graph (dict): The De Bruijn graph as adjacency lists.
    
    Returns:
        defaultdict: A mapping from nodes to their neighbors.
    """
    node_edge_map = defaultdict(list)
    for node, neighbors in graph.items():
        for neighbor in neighbors:
            node_edge_map[node].append(neighbor)
    return node_edge_map

def eulerian_trail(node_edge_map, start):
    """
    Finds an Eulerian trail in the De Bruijn graph.
    
    Args:
        node_edge_map (dict): The graph's node-to-edge mapping.
        start (str): The starting node for the trail.
    
    Returns:
        list: A list of nodes representing the Eulerian trail.
    """
    trail = [start]

    while node_edge_map:
        current_trail = []
        current = start
        while current in node_edge_map:
            next_node = node_edge_map[current].pop()
            if not node_edge_map[current]:
                del node_edge_map[current]
            current_trail.append(next_node)
            current = next_node
        trail.extend(current_trail)
        start = next((n for n in trail if n in node_edge_map), None)
        if not start:
            break
    return trail

def assemble(trail):
    """
    Assembles a SPSS sequence from an Eulerian trail.
    
    Args:
        trail (list): A list of nodes representing the Eulerian trail.
    
    Returns:
        str: The assembled SPSS sequence.
    """
    if not trail:
        return ""
    spss = trail[0]
    for node in trail[1:]:
        spss += node[-1]
    return spss

def verify_spss_kmers(spss, solid_kmers, k):
    """
    Verifies that all solid k-mers are present in the assembled SPSS.
    
    Args:
        spss (str): The assembled SPSS sequence.
        solid_kmers (set): A set of solid canonical k-mers.
        k (int): Length of the k-mers.
    """
    reconstructed_kmers = {spss[i:i+k] for i in range(len(spss) - k + 1)}
    missing_kmers = set(solid_kmers) - reconstructed_kmers
    if missing_kmers:
        print(f"Warning: The following k-mers are missing from the SPSS: {missing_kmers}")
    else:
        print("All solid k-mers are present in the SPSS.")

def verify_fm_index(fm_index, solid_kmers):
    """
    Verifies that all solid k-mers are correctly indexed in the FM-index.
    
    Args:
        fm_index (FmIndex): The FM-index object.
        solid_kmers (set): A set of solid canonical k-mers.
    """
    for kmer in solid_kmers:
        if fm_index.contains(kmer):
            print(f"K-mer {kmer} found in the FM-index.")
        else:
            print(f"Error: K-mer {kmer} not found in the FM-index.")

if __name__ == "__main__":
    solid_kmers = ["ATGG", "TGGC", "GGCC", "GCCA", "CCAA", "CAAC", "AACG"]

    # Build the De Bruijn graph
    graph = build_de_bruijn_graph(solid_kmers)
    print("De Bruijn Graph:")
    for node, neighbors in graph.items():
        print(f"{node} -> {neighbors}")

    # Find the Eulerian trail
    node_edge_map = make_node_edge_map(graph)
    start = next(iter(graph))  
    trail = eulerian_trail(node_edge_map, start)
    print("\nEulerian Trail:", trail)

    # Assemble the SPSS
    spss = assemble(trail)
    print("\nAssembled Sequence (SPSS):", spss)

    # Verify the SPSS
    verify_spss_kmers(spss, solid_kmers, 4)  

    # Build the FM-index
    fm_index = FmIndex(spss)
    print("\nFM-index built.")

    # Verify the FM-index
    verify_fm_index(fm_index, solid_kmers)

    # Save the FM-index
    output_file = "fm_index.dump"
    fm_index.save(output_file)
    print(f"\nFM-index saved to file: {output_file}")
