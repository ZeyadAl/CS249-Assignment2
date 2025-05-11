#!/usr/bin/env python3
"""
Minimal De Bruijn Graph assembler without classes.
Takes a FASTQ file, k-mer size, outputs GFA and FASTA contig.
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
    """Build de Bruijn graph: nodes are (k-1)-mers, edges represent k-mers."""
    graph = defaultdict(list)
    nodes = set()
    for seq in seqs:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            left = kmer[:-1]
            right = kmer[1:]
            graph[left].append(right)
            nodes.add(left)
            nodes.add(right)
    return graph, nodes


def write_gfa(graph, nodes, k, out_gfa):
    """Write graph in GFA format."""
    id_map = {node: f"S{idx}" for idx, node in enumerate(nodes)}
    overlap = k - 1
    with open(out_gfa, 'w') as f:
        # Header
        f.write("H\tVN:Z:1.0\n")
        # Segments
        for node in nodes:
            f.write(f"S\t{id_map[node]}\t{node}\n")
        # Links
        for left, rights in graph.items():
            for right in rights:
                f.write(f"L\t{id_map[left]}\t+\t{id_map[right]}\t+\t{overlap}M\n")


def main():
    parser = argparse.ArgumentParser(description="De Bruijn graph assembler")
    parser.add_argument('-i', '--input', required=True, help="Input FASTQ file")
    parser.add_argument('-k', type=int, required=True, help="k-mer size")
    parser.add_argument('-g', '--out_gfa', required=True, help="Output graph GFA")
    args = parser.parse_args()

    seqs = parse_fastq(args.input)
    graph, nodes = build_graph(seqs, args.k)
    write_gfa(graph, nodes, args.k, args.out_gfa)

if __name__ == '__main__':
    main()
