import logging as lg
from pathlib import Path

from .assessment import assess_canisregions
from .canisgenes import identify_canisgenes
from .canisregions import identify_canisregions
from .input import read_files

def run_padis(
        genomes_path: Path, genes_path: Path, pangenome_file: Path, 
        output_dir: Path, write_intervals: bool = False, debug: bool = False
        ) -> None:
    
    lg.basicConfig(
        level = lg.DEBUG if debug else lg.INFO,
        format = "[%(asctime)s] %(levelname)s: %(message)s",
        datefmt = "%d/%m %H:%M:%S")

    genomes_files = read_files(genomes_path)
    genes_files = read_files(genes_path)

    canisgenes_file = output_dir / "canisgenes.csv"
    intervals_file = output_dir / "intervals.csv" if write_intervals else None
    canisregions_file = output_dir / "canisregions.csv"
    canisregions_ass_file = output_dir / "canisresion_assessed.csv"

    lg.info("Starting phase 1: identification of candidate insertion sequence"
        "genes")
    identify_canisgenes(genes_files, pangenome_file, canisgenes_file,
        intervals_file)

    lg.info("Starting phase 2: extension to candidate insertion sequence "
        "regions")
    identify_canisregions(canisgenes_file, genomes_files, canisregions_file)

    lg.info("Starting phase 3: assessment of candidate insertion sequence"
        "regions")
    assess_canisregions(canisregions_file, genomes_files, canisregions_ass_file)
