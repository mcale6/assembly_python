from typing import List
from assembly_python.utils.graph import Graph
from assembly_python.utils.pathway import Pathway
import networkx as nx

def graph_is_isomorphic(g1: Graph, g2: Graph) -> bool:
    """Check if two graphs are isomorphic using NetworkX."""
    def node_match(n1, n2):
        return n1.get('color') == n2.get('color')
        
    def edge_match(e1, e2):
        return e1.get('color') == e2.get('color')
        
    return nx.is_isomorphic(
        g1.graph, 
        g2.graph,
        node_match=node_match,
        edge_match=edge_match
    )
    
def identify_matching_substructures(current_pathway: Pathway, 
                                  original_graph: Graph) -> List[Pathway]:
    """Find potential substructures that could be duplicated."""
    new_pathways = []
    remnant = current_pathway.remnant
    
    # Get connected components using NetworkX
    components = list(nx.connected_components(remnant.graph))
    
    for component in components:
        # Create subgraph for this component
        subgraph = remnant.graph.subgraph(component)
        
        # Find all possible subgraphs of appropriate size
        for size in range(2, len(subgraph.nodes())):
            # Get all connected subgraphs of this size
            for sub_nodes in nx.node_connected_component(subgraph, list(component)[0]):
                sub = subgraph.subgraph(sub_nodes)
                
                # Find matching subgraphs
                matches = find_matching_subgraphs(remnant.graph, sub)
                
                for match in matches:
                    # Create new pathway with this duplication
                    new_pathway = create_pathway_with_duplication(
                        current_pathway,
                        Graph(
                            vertices=list(sub.nodes()),
                            edges=list(sub.edges()),
                            vertex_colors=[sub.nodes[v].get('color') for v in sub.nodes()],
                            edge_colors=[sub.edges[e].get('color') for e in sub.edges()]
                        ),
                        Graph(
                            vertices=list(match.nodes()),
                            edges=list(match.edges()),
                            vertex_colors=[match.nodes[v].get('color') for v in match.nodes()],
                            edge_colors=[match.edges[e].get('color') for e in match.edges()]
                        )
                    )
                    new_pathways.append(new_pathway)
    
    return new_pathways

def find_matching_subgraphs(graph: nx.Graph, subgraph: nx.Graph) -> List[nx.Graph]:
    """Find all occurrences of subgraph in graph using NetworkX's isomorphism."""
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph, 
        subgraph,
        node_match=lambda n1, n2: n1.get('color') == n2.get('color'),
        edge_match=lambda e1, e2: e1.get('color') == e2.get('color')
    )
    
    return [graph.subgraph(match) for match in matcher.subgraph_isomorphisms_iter()]

def create_pathway_with_duplication(current_pathway: Pathway, 
                                  subgraph: Graph,
                                  matching_subgraph: Graph) -> Pathway:
    """Create new pathway by adding matching subgraphs."""
    new_pathway = current_pathway.copy()
    
    # Add subgraph to pathway
    new_pathway.pathway.append(subgraph)
    
    # Create duplicate record
    new_pathway.duplicates.append([
        tuple(e) for e in subgraph.edges + matching_subgraph.edges
    ])
    
    # Update remnant using NetworkX operations
    remnant_edges = set(new_pathway.remnant.edges) - set(matching_subgraph.edges)
    new_remnant = new_pathway.remnant.graph.edge_subgraph(remnant_edges)
    
    new_pathway.remnant = Graph(
        vertices=list(new_remnant.nodes()),
        edges=list(new_remnant.edges()),
        vertex_colors=[new_remnant.nodes[v].get('color') for v in new_remnant.nodes()],
        edge_colors=[new_remnant.edges[e].get('color') for e in new_remnant.edges()]
    )
    
    # Update atom equivalences using NetworkX node mapping
    vertex_map = {v: v for v in new_pathway.remnant.vertices}
    new_pathway.update_atom_equivalences(vertex_map)
    
    return new_pathway