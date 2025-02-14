from typing import Set, Optional, List
from dataclasses import dataclass
from enum import Enum
from rdkit import Chem
import networkx as nx

class EdgeType(Enum):
    SINGLE = "single"
    DOUBLE = "double" 
    TRIPLE = "triple"
    AROMATIC = "aromatic"
    
    @staticmethod
    def from_rdkit_bond(bond: Chem.Bond) -> 'EdgeType':
        """Convert RDKit bond type to our EdgeType."""
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

class Graph:
    """NetworkX graph wrapper with chemical functionality."""
    
    def __init__(self, graph: Optional[nx.Graph] = None, edge_types: Optional[List[EdgeType]] = None):
        """Initialize with optional existing NetworkX graph and edge types."""
        self.graph = graph if graph is not None else nx.Graph()
        self.edge_types = edge_types if edge_types is not None else [e for e in EdgeType]  # Default to all EdgeTypes
    
    @classmethod
    def from_mol_file(cls, mol_file: str, edge_types: Optional[List[EdgeType]] = None) -> 'Graph':
        """Create graph from MOL file with specified edge types."""
        mol = Chem.MolFromMolFile(mol_file, sanitize=True)
        if mol is None:
            raise ValueError(f"Failed to read MOL file: {mol_file}")
        return cls.from_rdkit_mol(mol, edge_types)
    
    @classmethod
    def from_smiles(cls, smiles: str, edge_types: Optional[List[EdgeType]] = None) -> 'Graph':
        """Create graph from SMILES string with specified edge types."""
        mol = Chem.MolFromSmiles(smiles, sanitize=True)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return cls.from_rdkit_mol(mol, edge_types)
    
    @classmethod 
    def from_inchi(cls, inchi: str, edge_types: Optional[List[EdgeType]] = None) -> 'Graph':
        """Create graph from InChI string with specified edge types."""
        mol = Chem.MolFromInchi(inchi, sanitize=True)
        if mol is None:
            raise ValueError(f"Invalid InChI: {inchi}")
        return cls.from_rdkit_mol(mol, edge_types)
    
    @classmethod
    def from_rdkit_mol(cls, mol: Chem.Mol, edge_types: Optional[List[EdgeType]] = None) -> 'Graph':
        """Convert RDKit Mol object to Graph with specified edge types."""
        graph = nx.Graph()
        edge_types = edge_types if edge_types is not None else [e for e in EdgeType]  # Default to all EdgeTypes
        
        # Add atoms with attributes
        for atom in mol.GetAtoms():
            graph.add_node(atom.GetIdx(), 
                         color=atom.GetSymbol(),
                         aromatic=atom.GetIsAromatic())
            
        # Identify aromatic rings
        aromatic_atoms = set(atom.GetIdx() for atom in mol.GetAtoms() 
                            if atom.GetIsAromatic())
        
        # Add bonds with attributes, considering edge types
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            
            # Determine bond type
            bond_type = None
            if (begin_idx in aromatic_atoms and 
                end_idx in aromatic_atoms and 
                bond.GetIsAromatic()):
                bond_type = EdgeType.AROMATIC
            else:
                rdkit_bond_type = bond.GetBondType()
                if rdkit_bond_type == Chem.BondType.SINGLE:
                    bond_type = EdgeType.SINGLE
                elif rdkit_bond_type == Chem.BondType.DOUBLE:
                    bond_type = EdgeType.DOUBLE
                elif rdkit_bond_type == Chem.BondType.TRIPLE:
                    bond_type = EdgeType.TRIPLE

            # Add edge only if its type is in the allowed edge_types
            if bond_type in edge_types:
                graph.add_edge(begin_idx, end_idx, color=bond_type.value)
                
        return cls(graph, edge_types)
    
    @classmethod
    def from_sdf(cls, sdf_file: str, edge_types: Optional[List[EdgeType]] = None) -> List['Graph']:
        """Create graphs from SDF file with specified edge types."""
        supplier = Chem.SDMolSupplier(sdf_file, sanitize=True)
        graphs = []
        for i, mol in enumerate(supplier):
            if mol is None:
                print(f"Warning: Failed to read molecule {i} from SDF")
                continue
            graphs.append(cls.from_rdkit_mol(mol, edge_types))
        return graphs
    
    @classmethod
    def from_file(cls, path: str, edge_types: Optional[List[EdgeType]] = None) -> 'Graph':
        """Read graph from text file with specified edge types."""
        graph = nx.Graph()
        with open(path) as f:
            next(f)  # Skip name line
            nodes = [int(x) for x in next(f).strip().split()]
            graph.add_nodes_from(nodes)
            
            edge_line = next(f).strip().split()
            edges = [(int(edge_line[i]), int(edge_line[i+1])) 
                    for i in range(0, len(edge_line), 2)]
            graph.add_edges_from(edges)
            
            vertex_colors = next(f).strip()
            if vertex_colors!= "!":
                for node, color in zip(nodes, vertex_colors.split()):
                    graph.nodes[node]['color'] = color
                    
            edge_colors = next(f).strip()
            if edge_colors!= "!":
                for edge, color in zip(edges, edge_colors.split()):
                    graph.edges[edge]['color'] = color
                    
        return cls(graph, edge_types)
    
    def copy(self) -> 'Graph':
        """Create deep copy of graph."""
        return Graph(self.graph.copy(), self.edge_types.copy())
    
    def __str__(self) -> str:
        """String representation."""
        output = []
        nodes = list(self.graph.nodes())
        edges = list(self.graph.edges())
        output.append(f"Nodes: {nodes}")
        output.append(f"Edges: {edges}")
        colors = nx.get_node_attributes(self.graph, 'color')
        if colors:
            output.append(f"Node colors: {[colors.get(n) for n in nodes]}")
        edge_colors = nx.get_edge_attributes(self.graph, 'color')
        if edge_colors:
            output.append(f"Edge colors: {[edge_colors.get(e) for e in edges]}")
        return "\n".join(output)
    
    def create_matcher(self, other: 'Graph') -> nx.isomorphism.GraphMatcher:
        """Create an isomorphism matcher between this graph and another."""
        return nx.isomorphism.GraphMatcher(
            self.graph,
            other.graph,
            node_match=lambda n1, n2: n1.get('color') == n2.get('color'),
            edge_match=lambda e1, e2: e1.get('color') == e2.get('color')
        )

    def is_isomorphic(self, other: 'Graph') -> bool:
        """Check if this graph is isomorphic to another graph."""
        matcher = self.create_matcher(other)
        return next(matcher.isomorphisms_iter(), None) is not None