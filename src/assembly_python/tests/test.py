from assembly_python.common.graph import Graph
from assembly_python.algorithms.assembly_search import assembly

def test_isomorphism():
    # Test 1: Simple path graphs with same colors
    g1 = Graph(
        vertices=[0, 1, 2],
        edges=[(0, 1), (1, 2)],
        vertex_colors=['C', 'C', 'C'],
        edge_colors=['single', 'double']
    )

    g2 = Graph(
        vertices=[5, 6, 7],
        edges=[(5, 6), (6, 7)],
        vertex_colors=['C', 'C', 'C'],
        edge_colors=['single', 'double']
    )

    print("Test 1 - Same colored paths:", g1.is_isomorphic(g2))  # True

    # Test 2: Different vertex colors
    g3 = Graph(
        vertices=[0, 1, 2],
        edges=[(0, 1), (1, 2)],
        vertex_colors=['C', 'O', 'C'],
        edge_colors=['single', 'double']
    )

    print("Test 2 - Different colors:", g1.is_isomorphic(g3))  # False

    # Test 3: Same colors and isomorphic topology (path vs. star are isomorphic for 3 nodes)
    g4 = Graph(
        vertices=[0, 1, 2],
        edges=[(0, 1), (0, 2)],
        vertex_colors=['C', 'C', 'C'],
        edge_colors=['single', 'double']
    )

    print("Test 3 - Different but isomorphic topology:", g1.is_isomorphic(g4))  # True

    # Test 4: Same structure and colors in different order
    g5 = Graph(
        vertices=[7, 8, 9],
        edges=[(8, 7), (9, 8)],
        vertex_colors=['C', 'C', 'C'],
        edge_colors=['double', 'single']
    )

    print("Test 4 - Reversed edges:", g1.is_isomorphic(g5))  # True

    # Test 5: Cyclic vs linear (non-isomorphic)
    g6 = Graph(
        vertices=[0, 1, 2],
        edges=[(0, 1), (1, 2), (2, 0)],
        vertex_colors=['C', 'C', 'C'],
        edge_colors=['single', 'double', 'single']
    )

    print("Test 5 - Cyclic vs linear:", g1.is_isomorphic(g6))  # False

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

def test_aspirin_assembly():
    """Test that aspirin assembly matches Go implementation."""
    # Create aspirin graph
    graph = MoleculeHandler.from_inchi(
        "InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)"
    )
    
    # Run assembly
    pathways = assembly(graph, variant="shortest")
    
    # Verify results match Go implementation
    assert len(pathways) == 1
    pathway = pathways[0]
    
    # Check assembly index
    assert assembly_index(pathway, graph) == 8
    
    # Check number of fragments
    assert len(pathway.pathway) == 3
    
    # Check duplicates found
    assert len(pathway.duplicates) == 3
    
    # Verify specific fragment structures
    fragments = pathway.pathway
    
    # Fragment 1: 3-carbon fragment
    assert len(fragments[0].vertices) == 3
    assert len(fragments[0].edges) == 2
    
    # Fragment 2: 4-atom fragment with oxygens
    assert len(fragments[1].vertices) == 4
    assert len(fragments[1].edges) == 3
    
    # Fragment 3: 3-carbon fragment
    assert len(fragments[2].vertices) == 3
    assert len(fragments[2].edges) == 2

def test_matching_go_cases():
    """Test additional cases from Go implementation."""
    test_cases = [
        ("testdata/graphs/square.txt", 2),
        ("testdata/graphs/triangle.txt", 2),
        ("testdata/graphs/hexagon.txt", 3),
    ]
    
    for file_path, expected_index in test_cases:
        graph = Graph.from_file(file_path)
        pathways = assembly(graph, variant="shortest")
        assert assembly_index(pathways[0], graph) == expected_index

def test_edge_cases():
    """Test edge cases from Go implementation."""
    # Empty graph
    graph = Graph(vertices=[], edges=[])
    assembly(graph)
        
    # Single edge
    graph = Graph(
        vertices=[0, 1],
        edges=[(0, 1)],
        vertex_colors=['C', 'C'],
        edge_colors=['single']
    )
    pathways = assembly(graph)
    assert len(pathways) == 1
    assert assembly_index(pathways[0], graph) == 0
    
    # Disconnected components
    graph = Graph.from_file("testdata/graphs/disconnected.txt")
    pathways = assembly(graph)
    assert len(pathways[0].pathway) > 1  # Should find fragments


if __name__ == "__main__":
    test_molecule_handler()
    test_isomorphism()
    test_edge_cases()
    test_matching_go_cases()