from typing import List
from assembly_python.common.graph import Graph
from assembly_python.common.pathway import Pathway
from assembly_python.algorithms.substructure_matching import identify_matching_substructures

def extend_pathway_shortest(current_pathway: Pathway,
                          best_pathways: List[Pathway],
                          original_graph: Graph) -> None:
    """Find shortest assembly pathway."""
    if best_pathways and best_pathways[0].assembly_index(original_graph) <= Pathway.best_assembly_index(original_graph, current_pathway):
        return
        
    new_pathways = identify_matching_substructures(current_pathway, original_graph)
    
    for new_pathway in new_pathways:
        if not best_pathways or new_pathway.assembly_index(original_graph) < best_pathways[0].assembly_index(original_graph):
            best_pathways.clear()
            best_pathways.append(new_pathway)
        extend_pathway_shortest(new_pathway, best_pathways, original_graph)

def extend_pathway_all_shortest(current_pathway: Pathway,
                              best_pathways: List[Pathway],
                              original_graph: Graph) -> None:
    """Find all shortest assembly pathways."""
    if best_pathways and best_pathways[0].assembly_index(original_graph) < current_pathway.assembly_index(original_graph):
        return
        
    new_pathways = identify_matching_substructures(current_pathway, original_graph)
    
    for new_pathway in new_pathways:
        new_index = new_pathway.assembly_index(original_graph)
        
        if not best_pathways:
            best_pathways.append(new_pathway)
        elif new_index < best_pathways[0].assembly_index(original_graph):
            best_pathways.clear()
            best_pathways.append(new_pathway)
        elif new_index == best_pathways[0].assembly_index(original_graph):
            best_pathways.append(new_pathway)
            
        extend_pathway_all_shortest(new_pathway, best_pathways, original_graph)

def extend_pathway_all(current_pathway: Pathway,
                      all_pathways: List[Pathway],
                      original_graph: Graph) -> None:
    """Find all possible assembly pathways."""
    new_pathways = identify_matching_substructures(current_pathway, original_graph)
    
    for new_pathway in new_pathways:
        all_pathways.append(new_pathway)
        extend_pathway_all(new_pathway, all_pathways, original_graph)