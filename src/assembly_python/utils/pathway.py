from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Set
from assembly_python.utils.graph import Graph
import networkx as nx

@dataclass
class Duplicates:
    """
    Represents duplicate fragments in the assembly pathway.
    left: Edges in the pathway fragment
    right: Corresponding edges in the remnant
    """
    left: List[Tuple[int, int]]
    right: List[Tuple[int, int]]
    atom_mapping: Dict[int, int] = field(default_factory=dict)  # Maps fragment atoms to remnant atoms

    def copy(self) -> 'Duplicates':
        """Create a deep copy of the Duplicates object."""
        return Duplicates(
            left=self.left.copy(),
            right=self.right.copy(),
            atom_mapping=self.atom_mapping.copy()
        )
    
@dataclass
class Pathway:
    """
    Represents an assembly pathway.
    pathway: List of graphs representing duplicated fragments
    remnant: Graph after removing duplicates 
    duplicates: List of duplicate edge pairs
    atom_equivalences: List of atom equivalence classes
    """
    pathway: List[Graph]
    remnant: Graph
    duplicates: List[List[int]]
    atom_equivalences: List[List[int]] = field(default_factory=list)

    def update_atom_equivalences(self, vertex_map: Dict[int, int], fragment_idx: int = -1) -> None:
        """
        Update atom equivalence classes based on a vertex mapping.
        
        Args:
            vertex_map: Dictionary mapping original vertex indices to new indices
            fragment_idx: Index of current fragment (-1 for remnant)
        """
        # Create new equivalence class if needed
        if fragment_idx == -1:
            self.atom_equivalences.append([])
            
        # Map vertices to their equivalence classes
        for orig_vertex, new_vertex in vertex_map.items():
            found = False
            for equiv_class in self.atom_equivalences:
                if orig_vertex in equiv_class:
                    equiv_class.append(new_vertex)
                    found = True
                    break
            if not found:
                self.atom_equivalences.append([orig_vertex, new_vertex])

    def format_atom_equivalences(self) -> str:
        """Format atom equivalences for output."""
        if not self.atom_equivalences:
            return "No atom equivalences"
        
        output = []
        for i, equiv_class in enumerate(self.atom_equivalences):
            output.append(f"Class {i+1}: {equiv_class}")
        return "\n".join(output)

    def copy(self) -> 'Pathway':
        """Create deep copy of pathway."""
        return Pathway(
            pathway=[g.copy() for g in self.pathway],
            remnant=self.remnant.copy(),
            duplicates=[d.copy() for d in self.duplicates],
            atom_equivalences=[ae.copy() for ae in self.atom_equivalences]
        )

    def __str__(self) -> str:
        """Format pathway for printing."""
        output = []
        output.append("Pathway Graphs:")
        for i, graph in enumerate(self.pathway):
            output.append(f"Graph {i+1}:")
            output.append(str(graph))
            output.append("------")
        
        output.append("\nRemnant Graph:")
        output.append(str(self.remnant))
        output.append("------")
        
        output.append("\nDuplicated Edges:")
        for i, dup in enumerate(self.duplicates):
            output.append(f"Duplicate {i+1}:")
            output.append(f"Left: {dup}")
        
        output.append("\nAtom Equivalences:")
        output.append(self.format_atom_equivalences())
            
        return "\n".join(output)

def format_pathway_output(pathway: Pathway, original_graph: Graph, verbose: bool = False) -> str:
    """Format pathway for human-readable output."""
    output = []
    
    # Add original graph info
    output.append("Original Graph:")
    output.append(str(original_graph))
    output.append("")
    
    # Add pathway fragments
    output.append("Pathway Fragments:")
    for i, fragment in enumerate(pathway.pathway):
        output.append(f"Fragment {i+1}:")
        output.append(str(fragment))
        output.append("")
    
    # Add remnant
    output.append("Remnant:")
    output.append(str(pathway.remnant))
    output.append("")
    
    # Add duplicates
    output.append("Duplicates:")
    for i, duplicate in enumerate(pathway.duplicates):
        output.append(f"Duplicate {i+1}: {duplicate}")
    output.append("")
    
    if verbose:
        # Add atom equivalences
        output.append("Atom Equivalences:")
        for i, equiv_class in enumerate(pathway.atom_equivalences):
            output.append(f"Class {i+1}: {equiv_class}")
        output.append("")
        
        # Add assembly index
        from assembly_python.algorithms.assembly_search import assembly_index
        output.append(f"Assembly Index: {assembly_index(pathway, original_graph)}")
        
    return "\n".join(output)