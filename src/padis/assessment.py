# terms/abbreviations: 
# - reg = region: DNA sequence around a gene to align 
# - ali = aligned section: part of a region that aligns to another region

from Bio import Align
from concurrent.futures import ProcessPoolExecutor
import logging as lg
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta, Sequence

from .input import read_canisgenes

#############
# top level #
#############

def assess_canisorthogroups(
        canisgenes_file: Path, assemblies_files: dict[str, Path], 
        canisorthogroups_file: Path, max_length: int, threads: int
        ) -> None:

    if canisorthogroups_file.exists():
        lg.info("Existing canisorthogroups file found - skipping phase 2")
        return() 
    
    lg.info("Reading candidate insertion sequence genes")
    canisgenes = read_canisgenes(canisgenes_file)

    lg.info("Indexing genome sequences (creates .fai files)")
    for p in assemblies_files.values():
        _ = Fasta(p)

    # # subset orthogroups for testing purposes
    # topn = canisgenes["orthogroup"].unique()[:20]
    # canisgenes = canisgenes[canisgenes["orthogroup"].isin(topn)]

    lg.info("Calculating basic orthogroup stats")
    canisogs1 = (
        canisgenes
        .groupby("orthogroup")
        .agg(
            genes = ("gene", "size"),
            genomes = ("genome", "nunique"),
            positions = ("position", "nunique")
        )
    )

    if threads == 1: 
        lg.info("Assessing flanking regions per orthogroup")
        canisogs2 = (
            canisgenes
            .groupby("orthogroup")
            .apply(
                lambda g: process_orthogroup(g, assemblies_files, max_length)
            )
        )
    else: 
        lg.info(
            f"Assessing flanking regions per orthogroup using {threads} "
            "threads")
        tasks = [(orthogroup, genes, assemblies_files, max_length) for
            orthogroup, genes in canisgenes.groupby("orthogroup", sort = False)]
        with ProcessPoolExecutor(max_workers = threads) as executor: 
            results = list(executor.map(_worker, tasks))
        canisogs2 = pd.DataFrame(results).set_index("orthogroup")

    lg.info("Merging orthogroup statistics")
    canisogs = (
        canisogs1
        .merge(canisogs2, on = "orthogroup")
        .reset_index()
    )

    lg.info("Writing candidate insertion sequence orthogroups")
    canisogs.to_csv(canisorthogroups_file, index = False)

###############
# lower level #
###############

def process_orthogroup(genes, assemblies_files, max_length):

    lg.debug(f"Processing orthogroup {genes.name}")

    result = pd.Series({
        "length": 0,
        "tir_score": 0,
        "fdr_score": 0,
        "problem": ""
    })

    # make copy to avoid modifying original gene table 
    genes = genes.copy()

    genes = genes[genes["position"].notna()]

    # open connections to assembly files
    # remark: indexing will not happen again since .fai files already exist
    genomes = {p.stem: Fasta(p) for p in assemblies_files.values()}

    gene1, reg1 = sequence_with_flanks(genomes, genes, max_length)
    genes = genes.loc[(genes["position"] != gene1.position)]
    gene2, reg2 = sequence_with_flanks(genomes, genes, max_length)

    if not reg1 or not reg2:
        result["problem"] = "All flanking regions too short"
        return(result)

    aligner = Align.PairwiseAligner(scoring = "blastn")
    aligner.mode = "local"
    strand = "+" if gene1.strand == gene2.strand else "-"
    alignments = aligner.align(reg1.seq, reg2.seq, strand = strand)
    alignment = alignments[0]
    alico = alignment.coordinates
    ali1_start, ali1_end = [reg1.start + c for c in sorted(alico[0, [0, -1]])]
    ali2_start, ali2_end = [reg2.start + c for c in sorted(alico[1, [0, -1]])]

    if (
        ali1_start > gene1.start
        or ali1_end < gene1.end
        or ali2_start > gene2.start
        or ali2_end < gene2.end
    ):
        result["problem"] = "Genes not fully inside the alignment"
        return(result)

    if (
        ali1_start == reg1.start
        or ali1_end == reg1.end
        or ali2_start == reg2.start
        or ali2_end == reg2.end
    ): 
        result["problem"] = "Alignment extends beyond flanks"
        return(result)

    # alico[0, 0] is always smaller than alico[0, -1]
    # (but if gene.strand == "-", alico[1, 0] is larger than alico[1, -1])
    ali1 = reg1[alico[0, 0]:alico[0, -1]]
    if gene1.strand == "-": ali1 = ali1.reverse.complement
    term_up = ali1[:30]
    term_down = ali1[-30:]
    term_down_rc = term_down.reverse.complement
    tir_alignments = aligner.align(term_up.seq, term_down_rc.seq)
    fdr_alignments = aligner.align(term_up.seq, term_down.seq)
    result["length"] = np.int64(len(ali1.seq))
    result["tir_score"] = np.int64(tir_alignments[0].score)
    result["fdr_score"] = np.int64(fdr_alignments[0].score)
    if len(ali1.seq) > max_length: 
        result["problem"] = "Aligned section too long"
    else:
        result["problem"] = ""
    return(result)

# helper for running process_orthogroup in parallel 
# needs to be top-level function
def _worker(task):
    orthogroup, genes, assemblies_files, max_length = task
    genes.name = orthogroup
    result = process_orthogroup(genes, assemblies_files, max_length)
    result.loc["orthogroup"] = orthogroup
    return(result)

def sequence_with_flanks(
        genomes: dict[Fasta], genes: pd.DataFrame, max_length: int
        ) -> tuple[tuple, Sequence]:
    """
    Sample a gene and extract its sequencing including flanking regions.

    If no complete flanks can be extracted or if many uncalled bases are found,
    a different gene will be sampled. If no genes meet these criteria, None will
    be returned as sequence. 
    
    :param genomes: Genome sequence (fna file) for every genome
    :param genes: Coordinates of all genes 
    :param flanksize: Size of the flanking region to extract
    """

    gene = None
    seq = None

    for gene in genes.itertuples():
        start = min(gene.end - max_length, gene.start)
        end = max(gene.start + max_length, gene.end)
        if start < 0: continue
        seq = genomes[gene.genome][gene.contig][start:end]
        gene_length = gene.end - gene.start
        if len(seq.seq) < (2 * max_length - gene_length): continue
        if seq.seq.count("N") > 10: continue
        break 

    return(gene, seq)
