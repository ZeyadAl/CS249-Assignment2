#!/usr/bin/env python3
import argparse
from collections import defaultdict, deque

def parse_fastq(fq_file):
    reads = {}
    with open(fq_file) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()
            f.readline()
            reads[header] = seq
    return reads

def overlap(a: str, b: str, min_length: int) -> int:
    start = 0
    while True:
        start = a.find(b[:min_length], start)
        if start == -1:
            return 0
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1

def build_overlap_graph(reads: dict[str, str], min_overlap: int) -> dict[str, list[tuple[str, int]]]:
    graph = defaultdict(list)
    for id_a, seq_a in reads.items():
        for id_b, seq_b in reads.items():
            if id_a == id_b:
                continue
            olen = overlap(seq_a, seq_b, min_overlap)
            if olen > 0:
                graph[id_a].append((id_b, olen))
    for node in graph:
        graph[node].sort(key=lambda x: x[1], reverse=True)
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

def assemble_contigs(paths: list[list[str]], reads: dict[str, str], graph: dict[str, list[tuple[str,int]]]) -> dict[str, str]:
    overlap_len = { (u,v):olen for u, edges in graph.items() for v, olen in edges }
    contigs = {}
    for idx, path in enumerate(paths, start=1):
        if not path:
            continue
        if len(path) == 1:
            contigs[f"contig_{idx}"] = reads[path[0]]
        else:
            seq = reads[path[0]]
            for i in range(len(path)-1):
                u, v = path[i], path[i+1]
                olen = overlap_len.get((u,v), 0)
                if olen == 0:
                    olen = overlap(reads[u], reads[v], 1)
                if olen > 0:
                    seq += reads[v][olen:]
                else:
                    seq += reads[v]
            contigs[f"contig_{idx}"] = seq
    return contigs

def write_fasta(contigs: dict[str,str], out_file: str):
    with open(out_file, 'w') as f:
        for cid, seq in contigs.items():
            f.write(f'>{cid}\n')
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def remove_transitive_edges(graph: dict[str, list[tuple[str,int]]]) -> dict[str, list[tuple[str,int]]]:
    graph = {u: list(edges) for u, edges in graph.items()}
    while True:
        neighbors = {u: {v for v, _ in edges} for u, edges in graph.items()}
        removed_one = False
        for u, edges in graph.items():
            for v, w in list(edges):
                for mid in neighbors[u]:
                    if mid != v and v in neighbors.get(mid, ()):
                        graph[u].remove((v, w))
                        removed_one = True
                        break
                if removed_one:
                    break
            if removed_one:
                break
        if not removed_one:
            break
    return graph

def find_max_weight_hamiltonian_path(graph: dict[str, list[tuple[str,int]]]) -> tuple[list[str], int]:
    """
    Held–Karp DP for exact max‐weight Hamiltonian path (feasible for n≲20).
    Returns (best_path, best_weight).
    """
    nodes = list(graph.keys())
    n = len(nodes)
    idx = {u:i for i,u in enumerate(nodes)}
    wmat = [[0]*n for _ in range(n)]
    for u, edges in graph.items():
        ui = idx[u]
        for v, weight in edges:
            wmat[ui][idx[v]] = weight

    Nmask = 1 << n
    dp = [[-float('inf')]*n for _ in range(Nmask)]
    parent = [[None]*n for _ in range(Nmask)]

    for i in range(n):
        dp[1<<i][i] = 0

    for mask in range(Nmask):
        for u in range(n):
            if not (mask & (1<<u)): 
                continue
            base = dp[mask][u]
            if base == -float('inf'):
                continue
            for v in range(n):
                if mask & (1<<v):
                    continue
                nm = mask | (1<<v)
                val = base + wmat[u][v]
                if val > dp[nm][v]:
                    dp[nm][v] = val
                    parent[nm][v] = u

    full = (1<<n) - 1
    best_w = max(dp[full])
    end = dp[full].index(best_w)

    path_idx = []
    mask, node = full, end
    while node is not None:
        path_idx.append(node)
        prev = parent[mask][node]
        mask ^= (1<<node)
        node = prev
    path_idx.reverse()
    return [nodes[i] for i in path_idx], best_w

def ont_find_hamiltonian_path(graph: dict[str, list[tuple[str,int]]], threshold=20) -> list[str]:
    """
    Wrapper: exact DP if |V|≤threshold, else greedy heuristic.
    """
    if len(graph) <= threshold:
        path, w = find_max_weight_hamiltonian_path(graph)
        return path

    best_path, best_w = [], -1
    for start in graph:
        visited = {start}
        path = [start]
        total = 0
        while True:
            choices = [(weight, nbr) for nbr, weight in graph[path[-1]] if nbr not in visited]
            if not choices:
                break
            weight, nbr = max(choices)
            path.append(nbr)
            visited.add(nbr)
            total += weight
        if total > best_w:
            best_w, best_path = total, path
    return best_path

def main():
    parser = argparse.ArgumentParser(description='OLC assembler optimized for error-free reads with max-weight Hamiltonian path')
    parser.add_argument('--ont', action='store_true',
                        help='Use ONT mode: remove transitive edges & exact/heuristic DP path')
    parser.add_argument('fastq', nargs='+', help='Input FASTQ files')
    parser.add_argument('-n', '--min_overlap', type=int, default=2, help='Minimum overlap length')
    parser.add_argument('-o', '--out', default='contigs.fasta', help='Output FASTA file')
    args = parser.parse_args()

    # 1) parse reads
    reads = {}
    for fq in args.fastq:
        reads.update(parse_fastq(fq))

    # 2) build and reduce overlap graph
    graph     = build_overlap_graph(reads, args.min_overlap)
    if args.ont:
        red_graph = remove_transitive_edges(graph)

        # 3) find best Hamiltonian path
        best_path = ont_find_hamiltonian_path(red_graph)

        # 4) assemble a single contig from that path
        contig_dict = assemble_contigs([best_path], reads, red_graph)

        # 5) write output
        write_fasta(contig_dict, args.out)
        length = len(next(iter(contig_dict.values())))
        print(f'Wrote 1 contig (length={length}) to {args.out}')
    else:
        #ham_path = find_hamiltonian_path(graph)
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
        
            
        # Assemble contigs from these paths
        contigs = assemble_contigs(paths, reads, graph)

    
        #contigs = assemble_contigs([ham_path], reads, graph)
    
        # Write contigs to FASTA file
        write_fasta(contigs, args.out)
    
        print(f'Wrote {len(contigs)} contigs to {args.out}')

if __name__ == '__main__':
    main()