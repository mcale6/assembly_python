from typing import List, Set
import networkx as nx  
from assembly_python.common.graph import Graph
from assembly_python.common.pathway import Pathway

def identify_matching_substructures(current_pathway: Pathway, original_graph: Graph) -> List[Pathway]:
    """Find all matching substructures and create new pathways."""
    new_pathways = []
    remnant = current_pathway.remnant
    
    # Get edge adjacencies like Go
    edge_adjacencies = {}
    edge_list = list(remnant.graph.edges())
    for i, (v1, v2) in enumerate(edge_list):
        adj = []
        for j, (u1, u2) in enumerate(edge_list):
            if i != j and (v1 in (u1, u2) or v2 in (u1, u2)):
                adj.append(j)
        edge_adjacencies[i] = sorted(adj)
    
    # For each starting edge
    for start_edge_idx in range(len(edge_list)):
        # Track edges in current subgraph
        sub_edges = [start_edge_idx]
        forbidden = {}  # edges we've tried
        
        while True:
            # Try to grow subgraph
            found_neighbor = False
            for edge_idx in sub_edges:
                for adj_edge_idx in edge_adjacencies[edge_idx]:
                    if not forbidden.get(adj_edge_idx, False) and adj_edge_idx not in sub_edges:
                        # Found new edge to add
                        sub_edges.append(adj_edge_idx)
                        found_neighbor = True
                        
                        # Get subgraph from these edges
                        sub_edge_pairs = [edge_list[i] for i in sub_edges]
                        nodes = set()
                        for e in sub_edge_pairs:
                            nodes.add(e[0])
                            nodes.add(e[1])
                        subgraph = Graph(remnant.graph.subgraph(nodes))
                        
                        # Look for matches
                        matcher = nx.algorithms.isomorphism.GraphMatcher(
                            remnant.graph, subgraph.graph,
                            node_match=lambda n1, n2: n1.get('color') == n2.get('color'),
                            edge_match=lambda e1, e2: e1.get('color') == e2.get('color')
                        )
                        
                        # Check each match
                        for mapping in matcher.subgraph_isomorphisms_iter():
                            matched_nodes = set(mapping.values())
                            if not matched_nodes.intersection(nodes):
                                new_pathway = current_pathway.create_with_duplication(
                                    subgraph,
                                    Graph(remnant.graph.subgraph(matched_nodes))
                                )
                                new_pathways.append(new_pathway)
                        break
                if found_neighbor:
                    break
            
            if not found_neighbor:
                # Backtrack
                if not sub_edges:
                    break
                last_edge = sub_edges.pop()
                forbidden[last_edge] = True
    
    return new_pathways
