from typing import List, Tuple, Dict, Set, Optional
from dataclasses import dataclass
from enum import Enum
from rdkit import Chem
import networkx as nx
import copy

class EdgeType(Enum):
    """
    Standardized edge types used throughout the application.
    Matches Go implementation's edge type strings.
    """
    SINGLE = "single"
    DOUBLE = "double"
    TRIPLE = "triple"
    AROMATIC = "aromatic"
    
    @staticmethod
    def from_rdkit_bond(bond: Chem.Bond) -> 'EdgeType':
        """Convert RDKit bond type to our standard edge type."""
        if bond.GetBondType() == Chem.BondType.SINGLE:
            return EdgeType.SINGLE
        elif bond.GetBondType() == Chem.BondType.DOUBLE:
            return EdgeType.DOUBLE
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            return EdgeType.TRIPLE
        elif bond.GetBondType() == Chem.BondType.AROMATIC:
            return EdgeType.AROMATIC
        else:
            raise ValueError(f"Unsupported bond type: {bond.GetBondType()}")
    
    @staticmethod
    def from_string(s: str) -> Optional['EdgeType']:
        """Convert string to edge type, case insensitive."""
        s = s.lower()
        for edge_type in EdgeType:
            if edge_type.value == s:
                return edge_type
        return None

class Graph:
    """Wrapper around NetworkX graph with additional chemical functionality"""
    def __init__(self, vertices=None, edges=None, vertex_colors=None, edge_colors=None):
        self.graph = nx.Graph()
        
        # Add vertices with colors
        if vertices:
            for i, v in enumerate(vertices):
                color = vertex_colors[i] if vertex_colors else None
                self.graph.add_node(v, color=color)
                
        # Add edges with colors
        if edges:
            for i, (v1, v2) in enumerate(edges):
                color = edge_colors[i] if edge_colors else None
                self.graph.add_edge(v1, v2, color=color)
                
    @property
    def vertices(self):
        return list(self.graph.nodes())
        
    @property 
    def edges(self):
        return list(self.graph.edges())
        
    @property
    def vertex_colors(self):
        return [self.graph.nodes[v].get('color') for v in self.vertices]
        
    @property
    def edge_colors(self):
        return [self.graph.edges[e].get('color') for e in self.edges]
        
    def copy(self):
        return Graph(
            vertices=self.vertices,
            edges=self.edges,
            vertex_colors=self.vertex_colors,
            edge_colors=self.edge_colors
        )
        
    def break_graph_on_edges(self, edge_indices: List[int]) -> Tuple['Graph', 'Graph']:
        """Split graph into two based on edge selection using NetworkX."""
        selected_edges = [self.edges[i] for i in edge_indices]
        break_graph = self.graph.edge_subgraph(selected_edges)
        remnant_edges = [e for i, e in enumerate(self.edges) if i not in edge_indices]
        remnant_graph = self.graph.edge_subgraph(remnant_edges)
        
        return (
            Graph(
                vertices=list(break_graph.nodes()),
                edges=list(break_graph.edges()),
                vertex_colors=[break_graph.nodes[v].get('color') for v in break_graph.nodes()],
                edge_colors=[break_graph.edges[e].get('color') for e in break_graph.edges()]
            ),
            Graph(
                vertices=list(remnant_graph.nodes()),
                edges=list(remnant_graph.edges()),
                vertex_colors=[remnant_graph.nodes[v].get('color') for v in remnant_graph.nodes()],
                edge_colors=[remnant_graph.edges[e].get('color') for e in remnant_graph.edges()]
            )
        )

    def get_edge_type(self, edge_idx: int) -> Optional[EdgeType]:
        """Get standardized edge type for given edge index."""
        if not self.edge_colors or edge_idx >= len(self.edge_colors):
            return None
        return EdgeType.from_string(self.edge_colors[edge_idx])

    def set_edge_type(self, edge_idx: int, edge_type: EdgeType) -> None:
        """Set standardized edge type for given edge index."""
        if not self.edge_colors:
            self.edge_colors = [EdgeType.SINGLE.value] * len(self.edges)
        if edge_idx >= len(self.edges):
            raise ValueError(f"Edge index {edge_idx} out of range")
        self.edge_colors[edge_idx] = edge_type.value

    def copy_with_edges(self, edge_indices: List[int]) -> 'Graph':
        """Create new graph copying only specified edges while preserving types."""
        new_edges = []
        new_edge_colors = [] if self.edge_colors else None
        
        for idx in edge_indices:
            if idx >= len(self.edges):
                raise ValueError(f"Edge index {idx} out of range")
            new_edges.append(self.edges[idx])
            if self.edge_colors:
                new_edge_colors.append(self.edge_colors[idx])
                
        # Get required vertices
        vertices_needed = set()
        for e1, e2 in new_edges:
            vertices_needed.add(e1)
            vertices_needed.add(e2)
            
        new_vertices = [v for v in self.vertices if v in vertices_needed]
        new_vertex_colors = (
            [c for v, c in zip(self.vertices, self.vertex_colors) if v in vertices_needed]
            if self.vertex_colors else None
        )
        
        return Graph(
            vertices=new_vertices,
            edges=new_edges,
            vertex_colors=new_vertex_colors,
            edge_colors=new_edge_colors,
        )

    def edge_adjacencies(self) -> Dict[int, List[int]]:
        """Create mapping of edge indices to their adjacent edge indices using NetworkX."""
        return {i: sorted([j for j, (u1, u2) in enumerate(self.edges) 
                if i != j and (v1 in (u1, u2) or v2 in (u1, u2))])
                for i, (v1, v2) in enumerate(self.edges)}

    def get_adjacency_list(self) -> Dict[int, List[int]]:
        """Create adjacency list representation using NetworkX."""
        return dict(self.graph.adjacency())

    def connected_component_edges(self) -> List[List[int]]:
        """Find connected components based on edges using NetworkX."""
        components = list(nx.connected_components(self.graph))
        edge_components = []
        for component in components:
            subgraph = self.graph.subgraph(component)
            edge_indices = []
            for edge in subgraph.edges():
                edge_indices.append(self.edges.index(edge))
            edge_components.append(sorted(edge_indices))
        return sorted(edge_components)

    @classmethod
    def from_file(cls, file_path: str) -> 'Graph':
        """Read graph from custom text format."""
        with open(file_path) as f:
            lines = f.readlines()
            
        # Skip name line
        vertices = [int(x) for x in lines[1].strip().split()]
        edges = []
        for i in range(0, len(lines[2].strip().split()), 2):
            v1, v2 = map(int, lines[2].strip().split()[i:i+2])
            edges.append((v1, v2))
            
        vertex_colors = None if lines[3].strip() == "!" else lines[3].strip().split()
        edge_colors = None if lines[4].strip() == "!" else lines[4].strip().split()
        
        return cls(vertices=vertices, edges=edges, 
                  vertex_colors=vertex_colors, edge_colors=edge_colors)

    def __str__(self) -> str:
        """String representation of graph."""
        output = []
        output.append(f"Vertices: {self.vertices}")
        output.append(f"Edges: {self.edges}")
        if self.vertex_colors:
            output.append(f"Vertex colors: {self.vertex_colors}")
        if self.edge_colors:
            output.append(f"Edge colors: {self.edge_colors}")
        return "\n".join(output)