from typing import List, Optional
from rdkit import Chem
import networkx as nx
from assembly_python.utils.graph import Graph
from enum import Enum
from assembly_python.utils.graph import EdgeType

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

class MoleculeHandler:
    """
    Handles molecule input from various formats and converts to Graph structure.
    Supports: SMILES, InChI, MOL files, SDF files, PDB files.
    """

    @staticmethod
    def from_smiles(smiles: str, sanitize: bool = True) -> Graph:
        """Convert SMILES string to Graph."""
        mol = Chem.MolFromSmiles(smiles, sanitize=sanitize)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {smiles}")
        return MoleculeHandler._mol_to_graph(mol)
    
    @staticmethod
    def from_inchi(inchi: str, sanitize: bool = True) -> Graph:
        """Convert InChI string to Graph."""
        mol = Chem.MolFromInchi(inchi, sanitize=sanitize)
        if mol is None:
            raise ValueError(f"Invalid InChI string: {inchi}")
        return MoleculeHandler._mol_to_graph(mol)
    
    @staticmethod
    def from_mol_file(file_path: str, sanitize: bool = True) -> Graph:
        """Read mol file and convert to Graph."""
        mol = Chem.MolFromMolFile(file_path, sanitize=sanitize)
        if mol is None:
            raise ValueError(f"Failed to read mol file: {file_path}")
        return MoleculeHandler._mol_to_graph(mol)
    
    @staticmethod
    def from_mol_block(mol_block: str, sanitize: bool = True) -> Graph:
        """Convert MOL block string to Graph."""
        mol = Chem.MolFromMolBlock(mol_block, sanitize=sanitize)
        if mol is None:
            raise ValueError(f"Invalid MOL block")
        return MoleculeHandler._mol_to_graph(mol)
    
    @staticmethod
    def from_sdf_file(file_path: str, sanitize: bool = True) -> List[Graph]:
        """Read SDF file and convert all molecules to Graphs."""
        supplier = Chem.SDMolSupplier(file_path, sanitize=sanitize)
        graphs = []
        for i, mol in enumerate(supplier):
            if mol is None:
                print(f"Warning: Failed to read molecule {i} from SDF")
                continue
            graphs.append(MoleculeHandler._mol_to_graph(mol))
        return graphs
    
    @staticmethod
    def from_pdb_file(file_path: str, sanitize: bool = True) -> Graph:
        """Read PDB file and convert to Graph."""
        mol = Chem.MolFromPDBFile(file_path, sanitize=sanitize)
        if mol is None:
            raise ValueError(f"Failed to read PDB file: {file_path}")
        return MoleculeHandler._mol_to_graph(mol)

    @staticmethod
    def _mol_to_graph(mol: Chem.Mol) -> Graph:
        """Convert RDKit Mol object to Graph matching Go implementation."""
        vertices = list(range(mol.GetNumAtoms()))
        vertex_colors = [atom.GetSymbol() for atom in mol.GetAtoms()]
        
        edges = []
        edge_colors = []
        
        # Identify aromatic rings using RDKit's built-in aromaticity detection
        aromatic_atoms = set(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic())
        
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            edges.append((begin_idx, end_idx))
            
            # Match Go implementation's aromatic bond detection
            if (begin_idx in aromatic_atoms and 
                end_idx in aromatic_atoms and 
                bond.GetIsAromatic()):
                edge_colors.append(EdgeType.AROMATIC.value)
            else:
                edge_colors.append(EdgeType.from_rdkit_bond(bond).value)
                
        return Graph(
            vertices=vertices,
            edges=edges,
            vertex_colors=vertex_colors,
            edge_colors=edge_colors
        )
