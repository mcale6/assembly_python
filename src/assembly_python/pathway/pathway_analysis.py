import math
from typing import List
from assembly_python.pathway import Pathway
from assembly_python.core.graph import Graph

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
    """Calculate assembly index for a pathway."""
    index = len(original_graph.edges) - 1
    index -= pathway_steps_saved(pathway, True)
    return index

def max_steps_saved(pathway: Pathway) -> int:
    """Calculate maximum possible additional steps that could be saved."""
    connected_components = pathway.remnant.connected_component_edges()
    max_steps_saved = 0

    for component in connected_components:
        num_edges = len(component)
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
    best_steps = pathway_steps_saved(best_pathways[0], True)
    new_steps = pathway_steps_saved(new_pathway, True)

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
    