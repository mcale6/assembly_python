from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Set
from assembly_python.common.graph import Graph

@dataclass
class Duplicates:
    left: List[tuple]  # Edges in the fragment
    right: List[tuple]  # Corresponding edges in remnant
    atom_mapping: Dict[int, int]  # Maps fragment atoms to remnant atoms

@dataclass
class Pathway:
    pathway: List[Graph]  # List of duplicated fragments
    remnant: Graph       # Remaining structure
    duplicates: List[Duplicates]  # List of duplicate records
    atom_equivalences: List[Set[int]] = field(default_factory=list)  # Sets of equivalent atoms

    def copy(self) -> 'Pathway':
        """Create deep copy."""
        return Pathway(
            pathway=[Graph(g.graph.copy()) for g in self.pathway],
            remnant=Graph(self.remnant.graph.copy()),
            duplicates=[Duplicates(
                left=d.left.copy(),
                right=d.right.copy(),
                atom_mapping=d.atom_mapping.copy()
            ) for d in self.duplicates],
            atom_equivalences=[s.copy() for s in self.atom_equivalences]
        )

    def update_atom_equivalences(self, new_mapping: Dict[int, int]) -> None:
        """Update atom equivalences with new mapping."""
        for atom1, atom2 in new_mapping.items():
            found = False
            for equiv_class in self.atom_equivalences:
                if atom1 in equiv_class or atom2 in equiv_class:
                    equiv_class.add(atom1)
                    equiv_class.add(atom2)
                    found = True
                    break
            if not found:
                self.atom_equivalences.append({atom1, atom2})

    def create_with_duplication(self, subgraph: Graph, matching_subgraph: Graph) -> 'Pathway':
        """Create new pathway with duplication."""
        new_pathway = self.copy()
        new_pathway.pathway.append(subgraph)
        
        # Use the Graph's create_matcher method
        matcher = matching_subgraph.create_matcher(subgraph)
        mapping = next(matcher.isomorphisms_iter())
        
        # Create duplicate record
        duplicate = Duplicates(
            left=list(subgraph.graph.edges()),
            right=[tuple(sorted([mapping[n1], mapping[n2]]))
                  for n1, n2 in subgraph.graph.edges()],
            atom_mapping=mapping
        )
        new_pathway.duplicates.append(duplicate)
        
        # Update remnant
        new_remnant = new_pathway.remnant.graph.copy()
        new_remnant.remove_edges_from(matching_subgraph.graph.edges())
        new_pathway.remnant = Graph(new_remnant)
        
        # Update atom equivalences
        new_pathway.update_atom_equivalences(mapping)
        
        return new_pathway

    def assembly_index(self, original_graph: Graph) -> int:
        """Calculate assembly index."""
        remnant_edges = len(self.remnant.graph.edges())
        pathway_edges = sum(len(fragment.graph.edges()) for fragment in self.pathway)
        duplicated_edges = sum(len(dup.left) for dup in self.duplicates)
        return remnant_edges + pathway_edges - duplicated_edges - 1

    @staticmethod
    def best_assembly_index(graph: Graph, current_pathway: 'Pathway') -> int:
        """Calculate best possible assembly index."""
        return len(graph.graph.edges()) - len(current_pathway.remnant.graph.edges())

def format_pathway_output(pathway: Pathway, original_graph: Graph, verbose: bool = False) -> str:
    """Format pathway for output."""
    output = []
    
    # Original graph
    output.append("Original Graph:")
    output.append(str(original_graph))
    output.append("")
    
    # Fragments
    output.append("Pathway Fragments:")
    for i, fragment in enumerate(pathway.pathway):
        output.append(f"Fragment {i+1}:")
        output.append(str(fragment))
        output.append("")
    
    # Remnant
    output.append("Remnant:")
    output.append(str(pathway.remnant))
    output.append("")
    
    # Duplicates
    output.append("Duplicates:")
    for i, duplicate in enumerate(pathway.duplicates):
        output.append(f"Duplicate {i+1}:")
        output.append(f"  Fragment edges: {duplicate.left}")
        output.append(f"  Remnant edges: {duplicate.right}")
        output.append(f"  Atom mapping: {duplicate.atom_mapping}")
    output.append("")
    
    # Atom equivalences
    output.append("Atom Equivalences:")
    for i, equiv_class in enumerate(pathway.atom_equivalences):
        output.append(f"Class {i+1}: {sorted(equiv_class)}")
    output.append("")
    
    if verbose:
        output.append(f"Assembly Index: {pathway.assembly_index(original_graph)}")
    
    return "\n".join(output)