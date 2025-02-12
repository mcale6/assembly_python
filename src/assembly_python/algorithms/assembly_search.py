from typing import List
import math
import networkx as nx
from assembly_python.utils.graph import Graph
from assembly_python.utils.pathway import Pathway
from assembly_python.algorithms.substructure_matching import identify_matching_substructures

def pathway_steps_saved(pathway: Pathway, edge_mode: bool = True) -> int:
    """Calculate steps saved in a pathway through reuse of fragments."""
    steps_saved = 0
    for fragment in pathway.pathway:
        if edge_mode:
            steps_saved += len(fragment.edges) - 1
        else:
            steps_saved += len(fragment.vertices) - 1
    return steps_saved

def assembly_index(pathway: Pathway, original_graph: Graph) -> int:
    """Calculate assembly index matching Go implementation."""
    # Count bonds in remnant
    remnant_bonds = len(pathway.remnant.edges)
    
    # Count bonds in pathway fragments
    pathway_bonds = sum(len(fragment.edges) for fragment in pathway.pathway)
    
    # Calculate duplicated bonds
    duplicated_bonds = sum(len(dup) for dup in pathway.duplicates)
    
    return remnant_bonds + pathway_bonds - duplicated_bonds - 1

def max_steps_saved(pathway: Pathway) -> int:
    """Calculate maximum possible additional steps that could be saved."""
    # Use NetworkX to find connected components
    components = list(nx.connected_components(pathway.remnant.graph))
    max_steps_saved = 0

    for component in components:
        # Get subgraph for this component
        subgraph = pathway.remnant.graph.subgraph(component)
        num_edges = len(subgraph.edges())
        naive_ma = num_edges - 1
        best_ma = int(math.floor(math.log2(float(num_edges))))
        best_steps_saved = naive_ma - best_ma
        max_steps_saved += best_steps_saved

    return max_steps_saved

def best_assembly_index(g: Graph, pathway: Pathway) -> int:
    """Calculate best possible assembly index based on maximum steps saved."""
    return assembly_index(pathway, g) - max_steps_saved(pathway)

def best_pathway_update(best_pathways: List[Pathway], new_pathway: Pathway,
                       variant: str = "shortest") -> None:
    """Update best pathways list based on variant strategy."""
    best_steps = pathway_steps_saved(best_pathways[0])
    new_steps = pathway_steps_saved(new_pathway)

    if variant == "shortest":
        if new_steps > best_steps:
            best_pathways.clear()
            best_pathways.append(new_pathway.copy())
    elif variant == "all_shortest":
        if new_steps > best_steps:
            best_pathways.clear()
            best_pathways.append(new_pathway.copy())
        elif new_steps == best_steps:
            best_pathways.append(new_pathway.copy())
    elif variant == "all":
        best_pathways.append(new_pathway.copy())
    else:
        raise ValueError(f"Unknown variant: {variant}")
    

def extend_pathway_shortest(current_pathway: Pathway,
                          best_pathways: List[Pathway],
                          original_graph: Graph) -> None:
    """Extend pathway looking for shortest assembly."""
    # Check if current pathway can be better than best found
    if best_pathways and assembly_index(best_pathways[0], original_graph) <= assembly_index(current_pathway, original_graph):
        return
        
    # Find substructures that could be duplicated
    new_pathways = identify_matching_substructures(current_pathway, original_graph)
    
    # Process new pathways
    for new_pathway in new_pathways:
        new_index = assembly_index(new_pathway, original_graph)
        
        # If better than current best, replace best
        if not best_pathways or new_index < assembly_index(best_pathways[0], original_graph):
            best_pathways.clear()
            best_pathways.append(new_pathway)
            
        # Continue search with new pathway
        extend_pathway_shortest(new_pathway, best_pathways, original_graph)


def extend_pathway_all_shortest(current_pathway: Pathway,
                              best_pathways: List[Pathway],
                              original_graph: Graph) -> None:
    """Extend pathway finding all shortest assemblies."""
    if best_pathways and assembly_index(best_pathways[0], original_graph) < assembly_index(current_pathway, original_graph):
        return
        
    new_pathways = identify_matching_substructures(current_pathway, original_graph)
    
    for new_pathway in new_pathways:
        new_index = assembly_index(new_pathway, original_graph)
        
        if not best_pathways:
            best_pathways.append(new_pathway)
        elif new_index < assembly_index(best_pathways[0], original_graph):
            best_pathways.clear()
            best_pathways.append(new_pathway)
        elif new_index == assembly_index(best_pathways[0], original_graph):
            best_pathways.append(new_pathway)
            
        extend_pathway_all_shortest(new_pathway, best_pathways, original_graph)


def extend_pathway_all(current_pathway: Pathway,
                      all_pathways: List[Pathway],
                      original_graph: Graph) -> None:
    """Extend pathway finding all possible assemblies."""
    new_pathways = identify_matching_substructures(current_pathway, original_graph)
    
    for new_pathway in new_pathways:
        all_pathways.append(new_pathway)
        extend_pathway_all(new_pathway, all_pathways, original_graph)