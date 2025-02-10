from typing import List
from rdkit import Chem
from assembly_python.core.graph import Graph

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
        """Read MOL file and convert to Graph."""
        mol = Chem.SDMolSupplier(file_path, sanitize=sanitize)[0]
        if mol is None:
            raise ValueError(f"Failed to read MOL file: {file_path}")
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
        """Convert RDKit Mol object to our Graph structure."""
        # Get vertices (atoms)
        vertices = list(range(mol.GetNumAtoms()))
        
        # Get edges (bonds)
        edges = []
        edge_colors = []
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            edges.append((begin_idx, end_idx))
            
            # Convert bond type to string representation
            bond_type = str(bond.GetBondType())
            if bond_type == "SINGLE":
                edge_colors.append("single")
            elif bond_type == "DOUBLE":
                edge_colors.append("double")
            elif bond_type == "TRIPLE":
                edge_colors.append("triple")
            elif bond_type == "AROMATIC":
                edge_colors.append("aromatic")
            else:
                edge_colors.append("unknown")
        
        # Get vertex colors (atom types)
        vertex_colors = [atom.GetSymbol() for atom in mol.GetAtoms()]
        
        return Graph(
            vertices=vertices,
            edges=edges,
            vertex_colors=vertex_colors,
            edge_colors=edge_colors
        )
