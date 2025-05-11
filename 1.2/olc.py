#!/usr/bin/env python3
import argparse
from collections import defaultdict, Counter, deque
import pprint

def parse_fastq(fq_file):
    """Parse a FASTQ file and return a dict of header->sequence."""
    reads = {}
    with open(fq_file) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # plus line
            f.readline()  # quality line
            reads[header] = seq
    return reads

def overlap(a: str, b: str, min_length: int) -> int:
    """
    Return length of longest suffix of 'a' matching a prefix of 'b'
    that is at least 'min_length'. Returns 0 if no such overlap.
    """
    start = 0  # start index for searching in a
    while True:
        # find occurrence of b's first min_length bases in a
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        # if the rest of a from 'start' matches the beginning of b
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def build_overlap_graph(reads: dict[str, str], min_overlap: int) -> dict[str, list[tuple[str, int]]]:
    """
    Given a dict of reads and a minimum overlap length, 
    return an overlap graph:
      { read_id: [(other_read_id, overlap_len), ...], ... }
    
    For error-free reads, we want to find all valid overlaps.
    """
    graph = defaultdict(list)
    
    # For each pair of reads, find their overlap
    for id_a, seq_a in reads.items():
        for id_b, seq_b in reads.items():
            if id_a == id_b:
                continue
                
            olen = overlap(seq_a, seq_b, min_overlap)
            if olen > 0:
                graph[id_a].append((id_b, olen))
    
    # Sort edges by overlap length (descending)
    for node in graph:
        graph[node].sort(key=lambda x: x[1], reverse=True)
    
    # Ensure all reads appear in the graph
    for read_id in reads:
        if read_id not in graph:
            graph[read_id] = []
            
    return dict(graph)

def find_hamiltonian_path(graph: dict[str, list[tuple[str, int]]]) -> list[str]:
    """
    For error-free reads with good coverage, we should be able to find a single path
    that visits each node exactly once (Hamiltonian path).
    
    This is an approximation algorithm that builds the path greedily.
    """
    if not graph:
        return []
        
    # Find potential start nodes (nodes with no incoming edges or more outgoing than incoming)
    in_degree = defaultdict(int)
    for node, edges in graph.items():
        for target, _ in edges:
            in_degree[target] += 1
    
    # Identify nodes with no incoming edges (potential start points)
    start_candidates = []
    for node in graph:
        if in_degree[node] == 0:
            start_candidates.append(node)
        elif in_degree[node] < len(graph.get(node, [])):
            # More outgoing than incoming edges
            start_candidates.append(node)
    
    # If no clear starting points, use any node
    if not start_candidates:
        start_candidates = list(graph.keys())
    
    # Try each potential start node
    best_path = []
    best_path_length = 0
    
    for start_node in start_candidates:
        # Start a path from this node
        path = [start_node]
        visited = {start_node}
        
        # Keep extending the path as long as possible
        current = start_node
        while True:
            # Find unvisited neighbors sorted by overlap length
            next_nodes = [(next_node, olen) for next_node, olen in graph.get(current, []) 
                          if next_node not in visited]
            
            if not next_nodes:
                break
                
            # Choose the neighbor with highest overlap
            next_nodes.sort(key=lambda x: x[1], reverse=True)
            next_node = next_nodes[0][0]
            
            # Add to path
            path.append(next_node)
            visited.add(next_node)
            current = next_node
        
        # Check if this path is better than our best so far
        if len(path) > best_path_length:
            best_path = path
            best_path_length = len(path)
        
        # If we've visited all nodes, we're done
        if len(visited) == len(graph):
            return path
    
    # After trying all start nodes, we have our best approximation
    remaining_nodes = set(graph.keys()) - set(best_path)
    
    # If there are remaining nodes, try to insert them where they fit best
    while remaining_nodes:
        best_insertion = None
        best_insertion_score = -1
        node_to_insert = None
        
        for node in remaining_nodes:
            # Try inserting this node at every position in the path
            for i in range(len(best_path)):
                # Check if this is a valid insertion point
                prev_node = best_path[i-1] if i > 0 else None
                next_node = best_path[i] if i < len(best_path) else None
                
                # Calculate insertion score based on overlaps
                score = 0
                if prev_node:
                    # Find overlap from prev_node to node
                    for target, olen in graph.get(prev_node, []):
                        if target == node:
                            score += olen
                            break
                
                if next_node:
                    # Find overlap from node to next_node
                    for target, olen in graph.get(node, []):
                        if target == next_node:
                            score += olen
                            break
                
                if score > best_insertion_score:
                    best_insertion = i
                    best_insertion_score = score
                    node_to_insert = node
        
        # If we found a good insertion point
        if best_insertion is not None:
            best_path.insert(best_insertion, node_to_insert)
            remaining_nodes.remove(node_to_insert)
        else:
            # If no good insertion points, just append to the path
            best_path.extend(list(remaining_nodes))
            break
    
    return best_path

def connect_paths(paths: list[list[str]], 
                 graph: dict[str, list[tuple[str, int]]],
                 reads: dict[str, str],
                 min_overlap: int) -> list[list[str]]:
    """
    Try to connect multiple paths into fewer paths by finding connections
    between the end of one path and the start of another.
    """
    if len(paths) <= 1:
        return paths
        
    # Create a graph of path connections
    path_graph = {}
    for i, path_i in enumerate(paths):
        path_graph[i] = []
        end_node_i = path_i[-1]
        
        for j, path_j in enumerate(paths):
            if i == j:
                continue
                
            start_node_j = path_j[0]
            
            # Check if there's a direct connection in the original graph
            direct_connection = False
            for target, _ in graph.get(end_node_i, []):
                if target == start_node_j:
                    direct_connection = True
                    break
            
            # Or compute overlap between end of path_i and start of path_j
            if not direct_connection:
                end_seq = reads[end_node_i]
                start_seq = reads[start_node_j]
                olen = overlap(end_seq, start_seq, min_overlap)
                direct_connection = (olen > 0)
            
            if direct_connection:
                path_graph[i].append(j)
    
    # Find connected components in the path graph
    visited = set()
    merged_paths = []
    
    for i in range(len(paths)):
        if i in visited:
            continue
            
        # BFS to find all connected paths
        component = []
        queue = deque([i])
        visited.add(i)
        
        while queue:
            node = queue.popleft()
            component.append(node)
            
            for neighbor in path_graph.get(node, []):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        
        # Sort component to get a valid path order
        ordered_component = []
        remaining = set(component)
        
        # Start with a node that has no incoming edges in this component
        start_candidates = list(remaining)
        for node in component:
            for neighbor in path_graph.get(node, []):
                if neighbor in remaining and neighbor in start_candidates:
                    start_candidates.remove(neighbor)
        
        # If no clear start, use any node
        current = start_candidates[0] if start_candidates else component[0]
        ordered_component.append(current)
        remaining.remove(current)
        
        # Build the ordered component
        while remaining:
            found = False
            for neighbor in path_graph.get(current, []):
                if neighbor in remaining:
                    ordered_component.append(neighbor)
                    remaining.remove(neighbor)
                    current = neighbor
                    found = True
                    break
            
            if not found:
                # If no direct neighbor, add any remaining node
                next_node = next(iter(remaining))
                ordered_component.append(next_node)
                remaining.remove(next_node)
                current = next_node
        
        # Merge the paths in this component
        merged_path = []
        for path_idx in ordered_component:
            if not merged_path:
                merged_path = paths[path_idx].copy()
            else:
                # Connect to the next path, avoid duplicating nodes if there's overlap
                next_path = paths[path_idx]
                
                # Check if the last node of merged_path connects to first node of next_path
                last_node = merged_path[-1]
                first_node = next_path[0]
                
                # If there's a direct connection in the original graph
                connection_exists = False
                for target, _ in graph.get(last_node, []):
                    if target == first_node:
                        connection_exists = True
                        break
                
                if connection_exists:
                    # Append the next path (without duplicating the first node)
                    merged_path.extend(next_path[1:])
                else:
                    # Just append the whole path
                    merged_path.extend(next_path)
        
        merged_paths.append(merged_path)
    
    return merged_paths

def assemble_contigs(
    paths: list[list[str]],
    reads: dict[str, str],
    graph: dict[str, list[tuple[str, int]]]
) -> dict[str, str]:
    """
    Given contig paths and the original reads + overlap graph, build
    each contig sequence by merging along the overlaps.
    """
    # Build quick lookup for overlap lengths
    overlap_len = {}
    for u, edges in graph.items():
        for v, olen in edges:
            overlap_len[(u, v)] = olen

    contigs = {}
    for idx, path in enumerate(paths, start=1):
        if not path:  # Skip empty paths
            continue
            
        if len(path) == 1:
            # Path with single read
            contigs[f"contig_{idx}"] = reads[path[0]]
        else:
            # Path with multiple reads
            seq = reads[path[0]]
            
            # For each pair of consecutive reads
            for i in range(len(path) - 1):
                prev, curr = path[i], path[i + 1]
                
                # If we have a stored overlap length, use it
                olen = overlap_len.get((prev, curr), 0)
                
                # If no stored overlap, compute it
                if olen == 0:
                    olen = overlap(reads[prev], reads[curr], 1)
                
                # Append only the non-overlapping suffix
                if olen > 0:
                    seq += reads[curr][olen:]
                else:
                    # If no overlap found, just append the whole sequence
                    # (this shouldn't happen with error-free reads, but as a safety measure)
                    seq += reads[curr]
            
            contigs[f"contig_{idx}"] = seq
    
    return contigs

def write_fasta(contigs, out_file):
    """Write contigs to a FASTA file."""
    with open(out_file, 'w') as f:
        for cid, seq in contigs.items():
            f.write(f'>{cid}\n')
            if cid == "contig_1":
                print(f"forst contig length is {len(seq)}")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def write_gfa(reads: dict[str, str],
              graph: dict[str, list[tuple[str,int]]],
              filename: str) -> None:
    """
    Write reads + overlap graph out in GFA v1 format.
    """
    with open(filename, 'w') as gfa:
        # Header
        gfa.write("H\tVN:Z:1.0\n")
        # Segments
        for rid, seq in reads.items():
            gfa.write(f"S\t{rid}\t{seq}\n")
        # Links
        for src, out_edges in graph.items():
            for dst, olen in out_edges:
                # '+' orientation for both ends, CIGAR = "<overlap>M"
                gfa.write(f"L\t{src}\t+\t{dst}\t+\t{olen}M\n")

def main():
    parser = argparse.ArgumentParser(description='OLC assembler optimized for error-free reads')
    parser.add_argument('fastq', nargs='+', help='Input FASTQ files')
    parser.add_argument('-n', '--min_overlap', type=int, default=2, help='Minimum overlap length')
    parser.add_argument('-o', '--out', default='contigs.fasta', help='Output FASTA file')
    parser.add_argument('--gfa', help='Optional GFA output file')
    parser.add_argument('--debug', action='store_true', help='Print debug information')
    args = parser.parse_args()

    # Parse reads from all input FASTQ files
    reads = {}
    for fq in args.fastq:
        reads.update(parse_fastq(fq))
    
    if args.debug:
        print(f"Loaded {len(reads)} reads")
    
    # Build the overlap graph
    graph = build_overlap_graph(reads, args.min_overlap)
    
    if args.debug:
        print("Original overlap graph:")
        edge_count = sum(len(edges) for edges in graph.values())
        print(f"Graph has {len(graph)} nodes and {edge_count} edges")
    
    # For error-free reads with good coverage, we should be able to find a Hamiltonian path
    # that visits each node exactly once
    ham_path = find_hamiltonian_path(graph)
    
    if args.debug:
        print(f"Hamiltonian path found with {len(ham_path)} nodes")
        if len(ham_path) == len(reads):
            print("Complete path found! This should generate 1 contig.")
        else:
            print(f"Incomplete path ({len(ham_path)}/{len(reads)} nodes). May generate multiple contigs.")
    
    # If we couldn't find a complete Hamiltonian path, fall back to the connect_paths approach
    if len(ham_path) < len(reads):
        # First try to find paths using classic OLC approach
        paths = []
        visited = set()
        
        # Start from nodes with no incoming edges
        in_degree = defaultdict(int)
        for node, edges in graph.items():
            for target, _ in edges:
                in_degree[target] += 1
        
        start_nodes = [node for node in graph if in_degree[node] == 0]
        
        # If no clear starting points, use nodes with high out-degree
        if not start_nodes:
            nodes_by_outdegree = sorted(graph.keys(), 
                                       key=lambda n: len(graph[n]), 
                                       reverse=True)
            start_nodes = nodes_by_outdegree[:max(1, len(nodes_by_outdegree) // 10)]
        
        # Create initial paths
        for start in start_nodes:
            if start in visited:
                continue
                
            path = [start]
            visited.add(start)
            current = start
            
            # Extend the path
            while True:
                next_nodes = [(node, olen) for node, olen in graph.get(current, []) 
                              if node not in visited]
                              
                if not next_nodes:
                    break
                    
                # Choose the highest overlap
                next_nodes.sort(key=lambda x: x[1], reverse=True)
                next_node = next_nodes[0][0]
                
                path.append(next_node)
                visited.add(next_node)
                current = next_node
            
            paths.append(path)
        
        # Handle unvisited nodes
        for node in graph:
            if node not in visited:
                paths.append([node])
                visited.add(node)
        
        # Try to connect these paths
        paths = connect_paths(paths, graph, reads, args.min_overlap)
        
        if args.debug:
            print(f"Found {len(paths)} paths after connection")
            
        # Assemble contigs from these paths
        contigs = assemble_contigs(paths, reads, graph)
    else:
        # We found a complete Hamiltonian path - just one contig!
        contigs = assemble_contigs([ham_path], reads, graph)
    
    # Write contigs to FASTA file
    write_fasta(contigs, args.out)
    
    print(f'Wrote {len(contigs)} contigs to {args.out}')
    #print(f'contig length is {len(contigs.items()[0][1])}')
    
    # Optionally write GFA file
    if args.gfa:
        write_gfa(reads, graph, args.gfa)
        print(f'Wrote overlap graph to {args.gfa}')

if __name__ == '__main__':
    main()
