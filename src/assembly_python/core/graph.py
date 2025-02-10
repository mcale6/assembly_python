from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass
from collections import defaultdict

@dataclass
class Graph:
    vertices: List[int]
    edges: List[Tuple[int, int]]
    vertex_colors: Optional[List[str]] = None
    edge_colors: Optional[List[str]] = None

    def __post_init__(self):
        """Validate graph initialization."""
        if not all(isinstance(v, int) for v in self.vertices):
            raise ValueError("All vertices must be integers")
        
        if not all(isinstance(e, tuple) and len(e) == 2 for e in self.edges):
            raise ValueError("Edges must be tuples of length 2")
        
        for v1, v2 in self.edges:
            if v1 not in self.vertices or v2 not in self.vertices:
                raise ValueError(f"Edge {(v1, v2)} references non-existent vertex")

        if self.vertex_colors and len(self.vertex_colors) != len(self.vertices):
            raise ValueError("Number of vertex colors must match number of vertices")
        
        if self.edge_colors and len(self.edge_colors) != len(self.edges):
            raise ValueError("Number of edge colors must match number of edges")

    def edge_adjacencies(self) -> Dict[int, List[int]]:
        """Map which edges share vertices."""
        out_map = defaultdict(list)
        for i, e1 in enumerate(self.edges):
            for j, e2 in enumerate(self.edges):
                if i != j and (e1[0] in e2 or e1[1] in e2):
                    out_map[i].append(j)
        return dict(out_map)

    def vertex_colour_map(self) -> Dict[int, str]:
        """Map vertices to their colors."""
        if not self.vertex_colors:
            return {}
        return {v: c for v, c in zip(self.vertices, self.vertex_colors)}

    def edge_colour_map(self) -> Dict[Tuple[int, int], str]:
        """Map edges to their colors."""
        if not self.edge_colors:
            return {}
        return {e: c for e, c in zip(self.edges, self.edge_colors)}

    def connected_component_edges(self) -> List[List[int]]:
        """Find sets of edges corresponding to all connected components."""
        edge_adj = self.edge_adjacencies()
        edges_used = set()
        edge_sets = []
        
        for i in range(len(self.edges)):
            if i not in edges_used:
                component = self.connected_component(i, edge_adj)
                edge_sets.append(component)
                edges_used.update(component)
        
        return edge_sets

    @staticmethod
    def connected_component(edge: int, edge_adj: Dict[int, List[int]]) -> List[int]:
        """Return connected component containing given edge."""
        component = [edge]
        i = 0
        while i < len(component):
            for j in edge_adj.get(component[i], []):
                if j not in component:
                    component.append(j)
            i += 1
        return component

    def copy_graph_edge(self, edge_index: int, new_graph: 'Graph') -> None:
        """Copy edge and associated vertices to new graph."""
        if edge_index >= len(self.edges):
            raise ValueError(f"Edge index {edge_index} out of range")

        edge = self.edges[edge_index]
        new_graph.edges.append(edge)

        if self.edge_colors:
            if not new_graph.edge_colors:
                new_graph.edge_colors = []
            new_graph.edge_colors.append(self.edge_colors[edge_index])

        for v in edge:
            if v not in new_graph.vertices:
                new_graph.vertices.append(v)
                if self.vertex_colors:
                    if not new_graph.vertex_colors:
                        new_graph.vertex_colors = []
                    v_idx = self.vertices.index(v)
                    new_graph.vertex_colors.append(self.vertex_colors[v_idx])

    def break_graph_on_edges(self, edge_indices: List[int]) -> Tuple['Graph', 'Graph']:
        """Split graph into two based on edge selection."""
        break_graph = Graph([], [])
        remnant_graph = Graph([], [])

        for i, edge in enumerate(self.edges):
            if i in edge_indices:
                self.copy_graph_edge(i, break_graph)
            else:
                self.copy_graph_edge(i, remnant_graph)

        return break_graph, remnant_graph

    @staticmethod
    def recombine_graphs(graph_left: 'Graph', graph_right: 'Graph') -> Tuple['Graph', Dict[int, int]]:
        """Combine two graphs, relabeling vertices of right graph."""
        output_edges = []
        output_vertices = []
        output_edge_colors = []
        output_vertex_colors = []

        # Copy left graph elements
        output_edges.extend(graph_left.edges)
        output_vertices.extend(graph_left.vertices)
        if graph_left.edge_colors:
            output_edge_colors.extend(graph_left.edge_colors)
        if graph_left.vertex_colors:
            output_vertex_colors.extend(graph_left.vertex_colors)

        # Calculate next vertex number
        next_vertex = max(max(graph_left.vertices), max(graph_right.vertices)) + 1
        vertex_map = {}

        # Copy right graph vertices with relabeling
        for i, vertex in enumerate(graph_right.vertices):
            new_vertex = vertex if vertex not in output_vertices else next_vertex
            if new_vertex == next_vertex:
                next_vertex += 1
            output_vertices.append(new_vertex)
            vertex_map[vertex] = new_vertex
            if graph_right.vertex_colors:
                output_vertex_colors.append(graph_right.vertex_colors[i])

        # Copy right graph edges with relabeled vertices
        for i, (v1, v2) in enumerate(graph_right.edges):
            output_edges.append((vertex_map[v1], vertex_map[v2]))
            if graph_right.edge_colors:
                output_edge_colors.append(graph_right.edge_colors[i])

        output_graph = Graph(
            vertices=output_vertices,
            edges=output_edges,
            vertex_colors=output_vertex_colors if output_vertex_colors else None,
            edge_colors=output_edge_colors if output_edge_colors else None
        )

        return output_graph, vertex_map

    def copy(self) -> 'Graph':
        """Create deep copy of graph."""
        return Graph(
            vertices=self.vertices.copy(),
            edges=self.edges.copy(),
            vertex_colors=self.vertex_colors.copy() if self.vertex_colors else None,
            edge_colors=self.edge_colors.copy() if self.edge_colors else None
        )