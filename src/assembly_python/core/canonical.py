from typing import List, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict
from assembly_python.core.graph import Graph

@dataclass
class AutomorphismData:
    """Tracks automorphisms found during canonicalization."""
    automorphisms: List[Graph]
    automorphism_levels: List[List[int]]

def individualize(partition: List[List[int]], vertex: int) -> List[List[int]]:
    """
    Individualizes a vertex within a partition.
    Example: individualize 3 in [[1], [2, 3, 4], [5, 6]] -> [[1], [3], [2, 4], [5, 6]]
    """
    individualized = []
    for part in partition:
        if vertex in part:
            # Add individualized vertex in its own part
            individualized.append([vertex])
            # Add remaining vertices in order
            remainder = [v for v in part if v != vertex]
            if remainder:
                individualized.append(remainder)
        else:
            individualized.append(part)
    return individualized

def degree_in_part(graph: Graph, vertex: int, part: List[int]) -> int:
    """Returns degree of vertex with respect to a set of vertices."""
    degree = 0
    for v1, v2 in graph.edges:
        if (vertex == v1 and v2 in part) or (vertex == v2 and v1 in part):
            degree += 1
    return degree

def shatter(graph: Graph, part_left: List[int], part_right: List[int]) -> List[List[int]]:
    """
    Split part_left into parts ordered by degree into part_right.
    Example: For square graph with part_left=(1,2,3) and part_right=(4),
    where 2 has degree 0 to (4) and both 1,3 have degree 1 to (4),
    returns [[2], [1, 3]]
    """
    degree_map = defaultdict(list)
    
    for v in part_left:
        degree = degree_in_part(graph, v, part_right)
        degree_map[degree].append(v)
    
    # Sort by degree and return parts
    degrees = sorted(degree_map.keys())
    return [degree_map[d] for d in degrees]

def is_equitable(graph: Graph, partition: List[List[int]]) -> bool:
    """
    Returns true if partition is equitable (vertices in same part have same degree to other parts).
    """
    for part_left in partition:
        for part_right in partition:
            if len(shatter(graph, part_left, part_right)) > 1:
                return False
    return True

def equitable_refinement(graph: Graph, partition: List[List[int]]) -> List[List[int]]:
    """
    Returns a more refined partition that is equitable.
    Continues splitting parts until no more splits are possible.
    """
    refined = [part.copy() for part in partition]
    
    while True:
        done = True
        for i, part in enumerate(refined):
            for j, target in enumerate(refined):
                shattering = shatter(graph, part, target)
                if len(shattering) > 1:
                    # Replace current part with its shattering
                    refined = refined[:i] + shattering + refined[i+1:]
                    done = False
                    break
            if not done:
                break
        if done:
            break
            
    return refined

def coarsest_equitable_colorings(graph: Graph, partition: List[List[int]]) -> Tuple[List[List[List[int]]], List[int]]:
    """
    Returns set of coarsest equitable colorings and individualized vertices.
    """
    if not is_equitable(graph, partition):
        return [equitable_refinement(graph, partition)], []
    
    output_colorings = []
    individualized_vertices = []
    
    # Find first non-trivial part (length > 1)
    for part in partition:
        if len(part) > 1:
            # Create colorings by individualizing each vertex
            for v in part:
                individualized = individualize(partition, v)
                refined = equitable_refinement(graph, individualized)
                output_colorings.append(refined)
                individualized_vertices.append(v)
            return output_colorings, individualized_vertices
            
    return [], []

def is_discrete(partition: List[List[int]]) -> bool:
    """Returns true if each part contains exactly one vertex."""
    return all(len(part) == 1 for part in partition)

def search_tree(graph: Graph, initial_coloring: Optional[List[List[int]]] = None, 
                prune: bool = True) -> Graph:
    """
    Main canonical graph algorithm, returns canonical version of input graph.
    """
    if initial_coloring is None:
        # Start with all vertices in one color class
        initial_coloring = [graph.vertices.copy()]
        
    canonical_container = [graph]  # Use list to allow modification in inner function
    automorphisms = AutomorphismData(automorphisms=[], automorphism_levels=[])
    tree_level = [0]
    backtrack = [False]
    backtrack_level = [-1]
    
    def search_tree_inner(colouring: List[List[int]], v: List[int]):
        """Recursive inner part of search tree algorithm."""
        if prune and backtrack[0]:
            if len(tree_level) == backtrack_level[0]:
                backtrack[0] = False
                backtrack_level[0] = -1
            else:
                return

        if is_discrete(colouring):
            # Create canonical candidate from discrete coloring
            perm = [part[0] for part in colouring]
            canonical_candidate = permute_graph(graph, perm)
            
            # Check for automorphism
            for i, g in enumerate(automorphisms.automorphisms):
                if graphs_equal(canonical_candidate, g):
                    backtrack[0] = True
                    backtrack_level[0] = max_equal_level(
                        automorphisms.automorphism_levels[i], 
                        tree_level)
                    return
                    
            # Add to automorphisms
            automorphisms.automorphisms.append(canonical_candidate)
            automorphisms.automorphism_levels.append(tree_level.copy())
            
            # Update canonical form if better
            if len(canonical_container[0].vertices) == 0 or \
               graph_greater_than(canonical_candidate, canonical_container[0]):
                canonical_container[0] = canonical_candidate
                
        else:
            # Get refinements and continue search
            colorings, vertices = coarsest_equitable_colorings(graph, colouring)
            
            for i, new_coloring in enumerate(colorings):
                new_v = v.copy()
                if vertices:
                    new_v.append(vertices[i])
                    
                tree_level.append(i)
                search_tree_inner(new_coloring, new_v)
                tree_level.pop()
                
    search_tree_inner(initial_coloring, [])
    return canonical_container[0]

def permute_graph(graph: Graph, permutation: List[int]) -> Graph:
    """Creates new graph with relabeled vertices according to permutation."""
    # Create vertex mapping
    vertex_map = {old: new for old, new in zip(sorted(graph.vertices), permutation)}
    
    # Create new edges with mapped vertices
    new_edges = [(vertex_map[v1], vertex_map[v2]) for v1, v2 in graph.edges]
    
    # Copy colors if present
    new_vertex_colors = graph.vertex_colors.copy() if graph.vertex_colors else None
    new_edge_colors = graph.edge_colors.copy() if graph.edge_colors else None
    
    return Graph(
        vertices=permutation.copy(),
        edges=new_edges,
        vertex_colors=new_vertex_colors,
        edge_colors=new_edge_colors
    )

def graph_greater_than(g1: Graph, g2: Graph) -> bool:
    """Compare graphs lexicographically by edge lists."""
    edges1 = sorted((min(e), max(e)) for e in g1.edges)
    edges2 = sorted((min(e), max(e)) for e in g2.edges)
    return edges1 > edges2

def graphs_equal(g1: Graph, g2: Graph) -> bool:
    """Check if two graphs are exactly equal (same vertices and edges)."""
    if set(g1.vertices) != set(g2.vertices):
        return False
    edges1 = {tuple(sorted(e)) for e in g1.edges}
    edges2 = {tuple(sorted(e)) for e in g2.edges}
    return edges1 == edges2

def max_equal_level(list1: List[int], list2: List[int]) -> int:
    """Returns index of last matching element between lists."""
    for i, (a, b) in enumerate(zip(list1, list2)):
        if a != b:
            return i - 1
    return min(len(list1), len(list2)) - 1

def graph_is_isomorphic(g1: Graph, g2: Graph) -> bool:
    """
    Determine if two graphs are isomorphic.
    Two graphs are isomorphic if they have:
    1. Same number of vertices and edges
    2. Same vertex colors (possibly in different order)
    3. Same edge colors (possibly in different order)
    4. Same connectivity pattern
    """
    # Basic checks first
    if len(g1.vertices) != len(g2.vertices) or len(g1.edges) != len(g2.edges):
        return False
        
    # Check colors if present
    if bool(g1.vertex_colors) != bool(g2.vertex_colors) or \
       bool(g1.edge_colors) != bool(g2.edge_colors):
        return False
        
    if g1.vertex_colors and sorted(g1.vertex_colors) != sorted(g2.vertex_colors):
        return False
    
    if g1.edge_colors and sorted(g1.edge_colors) != sorted(g2.edge_colors):
        return False
    
    # Create comparable degree sequences
    def get_vertex_signature(g: Graph, vertex: int) -> tuple:
        """Get a comparable signature for a vertex including its neighborhood."""
        neighbors = []
        vertex_idx = g.vertices.index(vertex)
        vertex_color = g.vertex_colors[vertex_idx] if g.vertex_colors else None
        
        for i, (v1, v2) in enumerate(g.edges):
            if v1 == vertex:
                neighbor_idx = g.vertices.index(v2)
                neighbor_color = g.vertex_colors[neighbor_idx] if g.vertex_colors else None
                edge_color = g.edge_colors[i] if g.edge_colors else None
                neighbors.append((neighbor_color, edge_color))
            elif v2 == vertex:
                neighbor_idx = g.vertices.index(v1)
                neighbor_color = g.vertex_colors[neighbor_idx] if g.vertex_colors else None
                edge_color = g.edge_colors[i] if g.edge_colors else None
                neighbors.append((neighbor_color, edge_color))
                
        return (vertex_color, tuple(sorted(neighbors)))
    
    # Get sorted vertex signatures for both graphs
    signatures1 = sorted(get_vertex_signature(g1, v) for v in g1.vertices)
    signatures2 = sorted(get_vertex_signature(g2, v) for v in g2.vertices)
    
    if signatures1 != signatures2:
        return False
        
    # If we got here, check full connectivity pattern
    def get_edge_signature(g: Graph) -> set:
        """Get a comparable signature for the graph's edge connectivity."""
        edges = set()
        for i, (v1, v2) in enumerate(g.edges):
            # Get vertex colors
            v1_color = g.vertex_colors[g.vertices.index(v1)] if g.vertex_colors else None
            v2_color = g.vertex_colors[g.vertices.index(v2)] if g.vertex_colors else None
            edge_color = g.edge_colors[i] if g.edge_colors else None
            
            # Create a canonical edge representation
            edge = tuple(sorted([(v1_color, 0), (v2_color, 1)]))  # Use 0,1 to break symmetry
            edges.add((edge, edge_color))
            
        return edges
        
    return get_edge_signature(g1) == get_edge_signature(g2)