from typing import List
import math
from assembly_python.core.graph import Graph
from assembly_python.pathway import Pathway, pathway_steps_saved, assembly_index, best_assembly_index, Duplicates
from assembly_python.core.canonical import graph_is_isomorphic

def check_subgraph_matches(current_pathway: Pathway, best_pathways: List[Pathway],
                          original_graph: Graph, subgraph: Graph, remnant: Graph) -> bool:
    """Check for matching subgraphs while only tracking best pathway."""
    k = len(subgraph.edges)
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
                        
                        # Replace best pathway if this one is better
                        best_steps = pathway_steps_saved(best_pathways[0], True)
                        new_steps = pathway_steps_saved(new_pathway, True)
                        if new_steps > best_steps:
                            best_pathways[0] = new_pathway
                            
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
    """Extend pathway looking for shortest assembly."""
    # Check if this pathway can potentially be better than best found
    if assembly_index(best_pathways[0], original_graph) < best_assembly_index(original_graph, current_pathway):
        return

    # Calculate sizes to check for subgraphs
    sizes_to_check = int(math.floor(float(len(current_pathway.remnant.edges)) / 2))

    # Get edge adjacencies for remnant
    edge_adjacencies = current_pathway.remnant.edge_adjacencies()

    # Initialize tracking structures
    forbidden = {}
    forbidden_size = {}
    
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