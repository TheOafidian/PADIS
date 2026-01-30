import argparse
from importlib.metadata import metadata
from pathlib import Path
import sys

from .padis import run_padis

def main() -> None:

    meta = metadata("padis")
    args = parse_arguments(meta) 

    intro = f"""Hi, this is {meta["Name"]} version {meta["Version"]}\n"""
    print(intro)

    genomes_path = Path(args.genomes)
    genes_path = Path(args.genes)
    pangenome_file = Path(args.pangenome)
    output_dir = Path(args.outputdir)

    run_padis(
        genomes_path, genes_path, pangenome_file, output_dir, 
        args.write_intervals, args.debug)

def parse_arguments(meta) -> None:

    parser = argparse.ArgumentParser(
        prog = meta["Name"].lower(), description = meta["Summary"])

    parser.add_argument(
        "genomes", 
        help = "genome sequences in fasta format, either as a file with paths"
        "or as a directory")
    parser.add_argument(
        "genes", 
        help = "gene coordinates in gff format, either as a file with paths"
        "or as a directory")
    parser.add_argument(
        "pangenome", help = "file with pangenome in SCARAP format")
    parser.add_argument(
        "outputdir", help = "directory for output")
    parser.add_argument(
        "-v", "--version", action = "version", 
        version = f"""{meta["Name"]} version {meta["Version"]}""")
    parser.add_argument(
        "-i", "--write_intervals", action = "store_true",
        help = "write left and right position indicators of intervals")
    parser.add_argument(
        "-b", "--debug", action = "store_true", 
        help = "be extra verbose for debugging")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    args = parser.parse_args() 

    return(args)

if __name__ == "__main__":
    main() 
