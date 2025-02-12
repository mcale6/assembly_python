from typing import List
import argparse
import sys
from pathlib import Path
from assembly_python.utils.graph import Graph
from assembly_python.utils.pathway import Pathway, format_pathway_output
from assembly_python.utils.molecules import MoleculeHandler
from assembly_python.algorithms.assembly_search import (
    extend_pathway_shortest,
    extend_pathway_all_shortest,
    extend_pathway_all,
)

def process_input_file(file_path: str) -> Graph:
    """Process input file and return corresponding graph."""
    ext = Path(file_path).suffix.lower()[1:]  # Remove dot from extension
    
    handlers = {
        'mol': MoleculeHandler.from_mol_file,
        'sdf': MoleculeHandler.from_sdf_file,
        'txt': Graph.from_file  # Added support for text graph format
    }
    
    if ext not in handlers:
        raise ValueError(f"Unsupported file format: {ext}")
        
    graph = handlers[ext](file_path)
    
    # Handle SDF files which return list of graphs
    if isinstance(graph, list):
        if not graph:
            raise ValueError("No valid molecules found in input file")
        graph = graph[0]
        
    return graph

def assembly(graph: Graph, variant: str = "shortest") -> List[Pathway]:
    """Main assembly algorithm."""
    if not graph.vertices or not graph.edges:
        raise ValueError("Invalid graph: Must have vertices and edges")
        
    if variant not in {"shortest", "all_shortest", "all"}:
        raise ValueError(f"Invalid variant: {variant}")
        
    # Initialize with empty pathway but proper remnant
    init_pathway = Pathway(
        pathway=[],
        remnant=graph.copy(),
        duplicates=[],
        atom_equivalences=[]
    )
    
    best_pathways = [init_pathway]
    
    # Map variant to appropriate algorithm
    algorithm_map = {
        "shortest": extend_pathway_shortest,
        "all_shortest": extend_pathway_all_shortest,
        "all": extend_pathway_all
    }
    
    extend_pathway = algorithm_map[variant]
    extend_pathway(init_pathway, best_pathways, graph)
    
    return best_pathways

def main() -> None:
    """Main entry point for CLI."""
    parser = argparse.ArgumentParser(description='Calculate Assembly Index for molecules')
    parser.add_argument('input_file', help='Input file path')
    parser.add_argument('--variant', default='shortest', choices=['shortest', 'all_shortest', 'all'],  help='Algorithm variant (default: shortest)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Output detailed pathway information')
    
    args = parser.parse_args()
    
    try:
        # Process input and run assembly
        graph = process_input_file(args.input_file)
        pathways = assembly(graph, variant=args.variant)
        
        if not pathways:
            print("No valid pathways found", file=sys.stderr)
            sys.exit(1)
            
        # Use format_pathway_output from pathway.py
        print(format_pathway_output(pathways[0], graph, args.verbose))
        
    except FileNotFoundError:
        print(f"Error: File not found - {args.input_file}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: Unexpected error - {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()