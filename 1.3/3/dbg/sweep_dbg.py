#!/usr/bin/env python3
"""
Unitig-based De Bruijn Graph assembler without classes.
Takes a FASTQ file, k-mer size, filters by abundance, outputs GFA and FASTA conitigs.
"""
import argparse
from collections import defaultdict


def parse_fastq(filename):
    """Yield sequences from FASTQ file."""
    with open(filename) as f:
        while True:
            header = f.readline().rstrip()
            if not header:
                break
            seq = f.readline().rstrip()
            f.readline()  # plus line
            f.readline()  # quality line
            yield seq


def build_graph(seqs, k):
    """Build de Bruijn graph with unique edges from k-mers seen >= min_count times."""
    # Count k-mers
    counts = defaultdict(int)
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            counts[seq[i:i+k]] += 1

    # Build adjacency with unique edges
    graph = defaultdict(list)
    for kmer, c in counts.items():
        left, right = kmer[:-1], kmer[1:]
        # add one edge per unique k-mer
        if right not in graph[left]:
            graph[left].append(right)
    return graph


def write_gfa(graph, k, out_gfa):
    """Write graph in GFA format, nodes are k-1-mers, links are overlaps."""
    nodes = list(set(graph.keys()) | {v for vs in graph.values() for v in vs})
    id_map = {node: f"S{idx}" for idx, node in enumerate(nodes)}
    overlap = k - 1
    with open(out_gfa, 'w') as f:
        f.write("H\tVN:Z:1.0\n")
        for node in nodes:
            f.write(f"S\t{id_map[node]}\t{node}\n")
        for left, rights in graph.items():
            for right in rights:
                f.write(f"L\t{id_map[left]}\t+\t{id_map[right]}\t+\t{overlap}M\n")


def get_conitigs(graph):
    """Generate maximal non-branching paths (conitigs)."""
    # compute degrees
    in_deg = defaultdict(int)
    out_deg = defaultdict(int)
    for u, vs in graph.items():
        out_deg[u] += len(vs)
        for v in vs:
            in_deg[v] += 1

    visited = set()
    conitigs = []
    # for each edge, if it starts a non-1-in-1-out node, walk forward
    for u, vs in graph.items():
        for v in vs:
            if (u, v) in visited:
                continue
            if in_deg[u] != 1 or out_deg[u] != 1:
                path = [u, v]
                visited.add((u, v))
                # extend forward while linear
                curr = v
                while in_deg[curr] == 1 and out_deg[curr] == 1:
                    nxt = graph[curr][0]
                    path.append(nxt)
                    visited.add((curr, nxt))
                    curr = nxt
                conitigs.append(path)
    return conitigs


def path_to_seq(path):
    """Convert a path of k-1-mers into assembled sequence."""
    seq = path[0]
    for node in path[1:]:
        seq += node[-1]
    return seq


def write_conitigs(conitigs, out_fasta, k):
    """Write conitigs as FASTA records."""
    with open(out_fasta, 'a') as f:
        total = 0
        for idx, path in enumerate(conitigs, 1):
            seq = path_to_seq(path)
            total += len(seq)
            #f.write(f">conitig_{idx} length={len(seq)}\n")
            #f.write(seq + "\n")

        f.write(f"k= {k}, num_contigs={idx}, contigs_length={total}\n")


def main():
    parser = argparse.ArgumentParser(description="Unitig-based De Bruijn assembler")
    parser.add_argument('-i', '--input', required=True, help="Input FASTQ file")
    parser.add_argument('-k', type=int, required=True, help="k-mer size")
    parser.add_argument('-o', '--out_fasta', required=True, help="Output conitigs FASTA file")
    #parser.add_argument('-g', '--out_gfa', required=True, help="Output graph GFA file")
    args = parser.parse_args()

    seqs = list(parse_fastq(args.input))
    graph = build_graph(seqs, args.k)
    #write_gfa(graph, args.k, args.out_gfa)
    conitigs = get_conitigs(graph)
    write_conitigs(conitigs, args.out_fasta, args.k)

if __name__ == '__main__':
    main()

