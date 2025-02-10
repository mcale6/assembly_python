from dataclasses import dataclass
from typing import List, Dict, Tuple
from assembly_python.core.graph import Graph

@dataclass
class Duplicates:
    """
    Represents duplicate fragments in the assembly pathway.
    left: Edges in the pathway fragment
    right: Corresponding edges in the remnant
    """
    left: List[Tuple[int, int]]
    right: List[Tuple[int, int]]

    def copy(self) -> 'Duplicates':
        """Create a deep copy of the Duplicates object."""
        return Duplicates(
            left=self.left.copy(),
            right=self.right.copy()
        )
    
@dataclass
class Pathway:
    """
    Represents an assembly pathway.
    pathway: List of graphs representing duplicated fragments
    remnant: Graph after removing duplicates 
    duplicates: List of duplicate edge pairs
    atom_equivalents: Groups of equivalent atoms across fragments
    """
    pathway: List[Graph]
    remnant: Graph
    duplicates: List[Duplicates]
    atom_equivalents: List[List[int]]

    def update_atom_equivalents(self, vertex_map: Dict[int, int]) -> None:
        """
        Update atom equivalence classes based on a vertex mapping.
        
        Args:
            vertex_map: Dictionary mapping original vertex indices to new indices
        """
        # Create new equivalence classes for mapped vertices
        for old_vertex, new_vertex in vertex_map.items():
            # Check if either vertex is already in an equivalence class
            old_class = None
            new_class = None
            
            for i, equiv_class in enumerate(self.atom_equivalents):
                if old_vertex in equiv_class:
                    old_class = i
                if new_vertex in equiv_class:
                    new_class = i
            
            if old_class is None and new_class is None:
                # Neither vertex is in an equivalence class, create new one
                self.atom_equivalents.append([old_vertex, new_vertex])
            elif old_class is None:
                # Add old vertex to new vertex's class
                self.atom_equivalents[new_class].append(old_vertex)
            elif new_class is None:
                # Add new vertex to old vertex's class
                self.atom_equivalents[old_class].append(new_vertex)
            elif old_class != new_class:
                # Merge the two equivalence classes
                self.atom_equivalents[old_class].extend(self.atom_equivalents[new_class])
                del self.atom_equivalents[new_class]

    def copy(self) -> 'Pathway':
        return Pathway(
            pathway=[g.copy() for g in self.pathway],
            remnant=self.remnant.copy(),
            duplicates=[d.copy() for d in self.duplicates],
            atom_equivalents=[ae.copy() for ae in self.atom_equivalents]
        )

    def __str__(self) -> str:
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
            output.append(f"Left: {dup.left}")
            output.append(f"Right: {dup.right}")
        
        output.append("\nAtom Equivalents:")
        for equiv in self.atom_equivalents:
            output.append(str(equiv))
            
        return "\n".join(output)