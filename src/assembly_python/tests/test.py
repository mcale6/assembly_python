from assembly_python.utils.molecules import MoleculeHandler
from assembly_python.core.graph import Graph
from assembly_python.core.canonical import graph_is_isomorphic

def test_isomorphism():
    """Test graph isomorphism with simple examples."""
    # Create two isomorphic graphs (square with different labelings)
    g1 = Graph(
        vertices=[1, 2, 3, 4],
        edges=[(1, 2), (2, 3), (3, 4), (4, 1)],
        vertex_colors=["C", "C", "C", "C"],
        edge_colors=["single"] * 4
    )
    
    g2 = Graph(
        vertices=[2, 3, 4, 1],
        edges=[(2, 3), (3, 4), (4, 1), (1, 2)],
        vertex_colors=["C", "C", "C", "C"],
        edge_colors=["single"] * 4
    )
    
    # Create non-isomorphic graph (different structure)
    g3 = Graph(
        vertices=[1, 2, 3, 4],
        edges=[(1, 2), (2, 3), (3, 4)],  # Missing one edge
        vertex_colors=["C", "C", "C", "C"],
        edge_colors=["single"] * 3
    )
    
    # Test isomorphism
    print("Testing graph isomorphism:")
    print(f"g1 and g2 isomorphic: {graph_is_isomorphic(g1, g2)}")  # Should be True
    print(f"g1 and g3 isomorphic: {graph_is_isomorphic(g1, g3)}")  # Should be False
    
    # Test with colored graphs
    g4 = Graph(
        vertices=[1, 2, 3, 4],
        edges=[(1, 2), (2, 3), (3, 4), (4, 1)],
        vertex_colors=["C", "N", "C", "N"],  # Alternating C/N colors
        edge_colors=["single"] * 4
    )
    
    g5 = Graph(
        vertices=[2, 3, 4, 1],
        edges=[(2, 3), (3, 4), (4, 1), (1, 2)],
        vertex_colors=["N", "C", "N", "C"],  # Same pattern, different ordering
        edge_colors=["single"] * 4
    )
    
    print(f"g4 and g5 isomorphic: {graph_is_isomorphic(g4, g5)}")  # Should be True

def test_molecule_handler():
    """Test the MoleculeHandler with various inputs."""
    print("Testing MoleculeHandler...")
    
    # Test SMILES input
    print("\nTesting SMILES input (ethanol):")
    try:
        graph = MoleculeHandler.from_smiles("CCO")
        print(f"Vertices: {graph.vertices}")
        print(f"Edges: {graph.edges}")
        print(f"Vertex colors: {graph.vertex_colors}")
        print(f"Edge colors: {graph.edge_colors}")
    except ValueError as e:
        print(f"Error: {e}")
    
    # Test InChI input
    print("\nTesting InChI input (aspirin):")
    try:
        graph = MoleculeHandler.from_inchi("InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)")
        print(f"Vertices: {graph.vertices}")
        print(f"Edges: {graph.edges}")
        print(f"Vertex colors: {graph.vertex_colors}")
        print(f"Edge colors: {graph.edge_colors}")
    except ValueError as e:
        print(f"Error: {e}")
    
    # Test MOL block input (if you have a MOL block string)
    print("\nTesting MOL block input (if available):")
    mol_block = """BENZOIC ACID, 2-(ACETYLOXY)-, ID: C50782
  NIST    15052607512D 1   1.00000     0.00000      
Copyright by the U.S. Sec. Commerce on behalf of U.S.A. All rights reserved.
 13 13  0     0  0              1 V2000
    1.7434    1.4944    0.0000 C   0  0  0  0  0  0           0  0  0
    1.7434    0.4981    0.0000 C   0  0  0  0  0  0           0  0  0
    0.8966    1.9925    0.0000 C   0  0  0  0  0  0           0  0  0
    2.5902    1.9925    0.0000 O   0  0  0  0  0  0           0  0  0
    0.8966    0.0000    0.0000 C   0  0  0  0  0  0           0  0  0
    0.0000    1.4944    0.0000 C   0  0  0  0  0  0           0  0  0
    0.8966    2.9887    0.0000 C   0  0  0  0  0  0           0  0  0
    0.0000    0.4981    0.0000 C   0  0  0  0  0  0           0  0  0
    0.0000    3.4869    0.0000 O   0  0  0  0  0  0           0  0  0
    1.7434    3.4869    0.0000 O   0  0  0  0  0  0           0  0  0
    3.4869    1.4944    0.0000 C   0  0  0  0  0  0           0  0  0
    4.3337    1.9925    0.0000 O   0  0  0  0  0  0           0  0  0
    3.4869    0.4981    0.0000 C   0  0  0  0  0  0           0  0  0
  1  2  2  0     0  0
  3  1  1  0     0  0
  1  4  1  0     0  0
  2  5  1  0     0  0
  6  3  2  0     0  0
  3  7  1  0     0  0
  4 11  1  0     0  0
  5  8  2  0     0  0
  8  6  1  0     0  0
  7  9  2  0     0  0
  7 10  1  0     0  0
 11 12  2  0     0  0
 11 13  1  0     0  0
M  END"""
    try:
        graph = MoleculeHandler.from_mol_block(mol_block)
        print(f"Vertices: {graph.vertices}")
        print(f"Edges: {graph.edges}")
        print(f"Vertex colors: {graph.vertex_colors}")
        print(f"Edge colors: {graph.edge_colors}")
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_molecule_handler()
    test_isomorphism()