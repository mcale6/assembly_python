from typing import List
import math
from assembly_python.core.graph import Graph
from assembly_python.pathway import Pathway, assembly_index, best_assembly_index, best_pathway_update, Duplicates
from assembly_python.core.canonical import graph_is_isomorphic
from assembly_python.core.canonical import graph_is_isomorphic

def check_subgraph_matches(current_pathway: Pathway, best_pathways: List[Pathway],
                          original_graph: Graph, subgraph: Graph, remnant: Graph) -> bool:
    """Check for matching subgraphs and track all shortest pathways."""
    k = len(subgraph.edges)  # Size of subgraph to match
    match = False
    
    # Get edge adjacencies for remnant
    edge_adjacencies = remnant.edge_adjacencies()
    
    # Initialize tracking structures
    forbidden = {}
    forbidden_size = {}
    sub = []
    
    # Process all edges in remnant
    for i in range(len(remnant.edges)):
        sub = [i]
        
        while True:
            # Find non-forbidden neighbor
            neighbor_found = False
            for e in sub:
                for adj in edge_adjacencies.get(e, []):
                    if not forbidden.get(adj, False) and adj not in sub:
                        if len(sub) <= k:
                            sub.append(adj)
                            neighbor_found = True
                            break
                if neighbor_found:
                    break
                    
            if neighbor_found:
                # Check subgraph match if size matches
                if len(sub) == k:
                    possible_duplicate, new_remnant = remnant.break_graph_on_edges(sub)
                    
                    # Check if subgraphs are isomorphic
                    if graph_is_isomorphic(subgraph, possible_duplicate):
                        match = True
                        
                        # Create new pathway with this match
                        new_pathway = current_pathway.copy()
                        new_pathway.pathway.append(subgraph.copy())
                        
                        # Add duplicate record using Duplicates class
                        new_duplicate = Duplicates(
                            left=subgraph.edges.copy(),
                            right=possible_duplicate.edges.copy()
                        )
                        new_pathway.duplicates.append(new_duplicate)
                        
                        # Update remnant and atom equivalents
                        new_graph, vertex_map = Graph.recombine_graphs(
                            new_remnant, possible_duplicate)
                        new_pathway.remnant = new_graph
                        new_pathway.update_atom_equivalents(vertex_map)
                        
                        # Update best pathways list - keeps all equally good pathways
                        best_pathway_update(best_pathways, new_pathway, "all_shortest")
                            
                        # Continue search with new pathway
                        extend_pathway(new_pathway, best_pathways, original_graph)
                continue
                
            # Backtracking
            if len(sub) == 0:
                break
                
            # Update forbidden status based on size comparisons
            this_forbid_size = len(sub)
            this_forbid = sub[-1]
            sub = sub[:-1]
            
            for k, v in forbidden_size.items():
                if v > this_forbid_size:
                    forbidden[k] = False
            forbidden[this_forbid] = True
            forbidden_size[this_forbid] = this_forbid_size

    return match

def extend_pathway(current_pathway: Pathway, best_pathways: List[Pathway],
                  original_graph: Graph) -> None:
    """Extend pathway finding all shortest assemblies."""
    # Check if this pathway can potentially match shortest found
    if (len(best_pathways) > 0 and 
        assembly_index(best_pathways[0], original_graph) < best_assembly_index(original_graph, current_pathway)):
        return

    # Calculate sizes to check for subgraphs
    sizes_to_check = int(math.floor(float(len(current_pathway.remnant.edges)) / 2))

    # Get edge adjacencies for remnant
    edge_adjacencies = current_pathway.remnant.edge_adjacencies()

    # Initialize tracking structures
    forbidden = {}
    forbidden_size = {}
    sub = []
    
    # Process all edges in remnant
    for i in range(len(current_pathway.remnant.edges)):
        sub = [i]
        
        while True:
            # Find non-forbidden neighbor
            neighbor_found = False
            for e in sub:
                for adj in edge_adjacencies.get(e, []):
                    if not forbidden.get(adj, False) and adj not in sub:
                        if len(sub) <= sizes_to_check:
                            sub.append(adj)
                            neighbor_found = True
                            break
                if neighbor_found:
                    break
                    
            if neighbor_found:
                # Extract and check subgraph
                subgraph, remnant = current_pathway.remnant.break_graph_on_edges(sub)
                if len(sub) > 1:
                    match = check_subgraph_matches(current_pathway, best_pathways,
                                                 original_graph, subgraph, remnant)
                    if match:
                        continue
                continue
                
            # Backtracking
            if len(sub) == 0:
                break
                
            # Update forbidden information
            this_forbid_size = len(sub)
            this_forbid = sub[-1]
            sub = sub[:-1]
            
            # Update forbidden lists
            for k, v in forbidden_size.items():
                if v > this_forbid_size:
                    forbidden[k] = False
            forbidden[this_forbid] = True 
            forbidden_size[this_forbid] = this_forbid_size