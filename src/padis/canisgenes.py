import logging as lg
import pandas as pd
from pathlib import Path

from .input import read_genes, read_pangenome

#############
# top level #
#############

def identify_canisgenes(
        genes_files: list[Path], pangenome_file: Path, canisgenes_file: Path,
        intervals_file: Path = None) -> None:
    """
    Identify candidate insertion sequence genes. 
    
    :param genes_files: Annotation files (gff). 
    :param pangenome_file: File with pangenome in SCARAP format (tsv file 
        without header with the columns gene, genome and orthogroup)
    :param canisgenes_file: Path to output file for candidate insertion sequence 
        genes. 
    :param intervals_file: Path to output file for intervals (left and right
        position indicators). 
    """
    
    lg.info("Reading pangenome")
    pan = read_pangenome(pangenome_file)
    genomes = pan["genome"].unique()
    orthogroups = pan["orthogroup"].unique()
    lg.info(
        f"Pangenome contains {len(pan)} genes, {len(genomes)} genomes "
        f"and {len(orthogroups)} orthogroups")

    lg.info("Identifying single-copy core orthogroups")
    core = determine_core(pan) 
    lg.info(f"Identified {len(core)} single-copy core orthogroups")
    
    lg.info("Initializing genomic positions")
    # "posind" stands for position indicator
    posinds = [o + "-" for o in core] + [o + "+" for o in core]
    positions = pd.DataFrame(
        {"posind": posinds, "position": range(len(posinds))})
    positions = positions.set_index("posind")

    lg.info("Extracting genome names from gene annotation paths")
    # why not create paths from genome names? 
    # --> paths may have different extensions and may be compressed
    genes_files_dict = {p.name[0:15]: p for p in genes_files}

    lg.info("Processing gene annotation files")
    pan["contig"] = pd.Series(dtype = "str")
    pan["strand"] = pd.Series(dtype = "str")
    # use pandas type "Int64" instead of "int" to be able to hold missing values
    pan["start"] = pd.Series(dtype = "Int64")
    pan["end"] = pd.Series(dtype = "Int64")
    pan["posind"] = pd.Series(dtype = 'str')
    # it is faster not to sort the index (pan = pan.sort_index())!
    # probable reason: for a unique index, pandas uses a hashmap for lookup
    # while for a sorted index it uses binary search
    pan = pan.set_index("gene")
    if intervals_file:
        with open(intervals_file, "w"):
            pass
    for genome in genomes: 
        lg.debug(f"Processing genome {genome}")
        genes = read_genes(genes_files_dict[genome])
        # the following does three things: 
        # - add gene coordinates to pangenome table
        # - merge positions of position indicators if they are adjacent
        # - add position indicators of accessory genes to pangenome table
        positions, pan = _process_genes(
            positions, pan, genes, core, intervals_file)
    positions = positions.reset_index() 
    pan = pan.reset_index() 
    pan = pan[pan["contig"].notnull()] # keep only accessory genes

    lg.info("Determining genomic positions of genes from position indicators")
    positions = dict(zip(positions["posind"], positions["position"]))
    pan["position"] = pan["posind"].map(positions).astype("Int64")
    pan = pan.drop("posind", axis = 1)
    lg.info(f"Identified {pan['position'].nunique()} unique genomic positions")
    lg.info(f"Assigned a position to {pan['position'].count().sum()} out of "
        f"{len(pan)} accessory genes")

    lg.info("Selecting candidate insertion sequence genes")
    canisgenes = pan.groupby("orthogroup").filter(
        lambda g: g["position"].nunique() > 1)
    lg.info(f"Selected {len(canisgenes)} candidate insertion sequence genes")

    lg.info("Writing candidate insertion sequence genes")
    canisgenes.to_csv(canisgenes_file, index = False)

###############
# lower level #
###############

def determine_core(pangenome: pd.DataFrame) -> set[str]:
    """
    Identify the single-copy core orthogroups in a pangenome. 
    
    :param pangenome: DataFrame with columns gene, genome and orthogroup. 
    :return: Set with single-copy core orthogroups. 
    """
    n_genomes = pangenome["genome"].nunique()
    core_orthogroups = pangenome.\
        value_counts(["genome", "orthogroup"]).\
        reset_index(name = "count").\
        query("count == 1").\
        drop("count", axis = 1).\
        value_counts("orthogroup").\
        reset_index(name = "count").\
        query(f"count == {n_genomes}")["orthogroup"].\
        unique()
    core_orthogroups = set(core_orthogroups)
    return(core_orthogroups)

def _process_genes(positions: pd.DataFrame, pan: pd.DataFrame, 
        genes: pd.DataFrame, core: set[str], intervals_file: Path = None
        ) -> (pd.DataFrame, pd.DataFrame):
    """
    Process gene annotation of a genome (helper of identify_canisgenes). 
    
    :param positions: DataFrame with posind and position, indexed on posind. 
    :param pan: DataFrame with gene, genome, orthogroup, contig, strand, start,
        end, and posind. Indexed on gene. 
    :param genes: DataFrame with columns of gff file. 
    :param core: Set with core orthogroups. 
    :param intervals_file: Path to intervals output file. 
    :return: Updated positions and pan DataFrames. 
    """
    
    contig = None
    interval = []
    left_posind = None
    right_posind = None

    # column indices for iloc (slightly faster than loc)
    genome_col = pan.columns.get_loc("genome")
    orthogroup_col = pan.columns.get_loc("orthogroup")
    contig_col = pan.columns.get_loc("contig")
    strand_col = pan.columns.get_loc("strand")
    start_col = pan.columns.get_loc("start")
    end_col = pan.columns.get_loc("end")
    
    # itertuples is faster than iterrows (tested)
    for row in genes.itertuples():

        # reconstruct gene id that prodigal uses in ffn/faa files but not
        # in gff file
        gene = row.seqid + "_" + row.attr.split(";")[0].split("_")[1]

        # row index for iloc (slightly faster than loc)
        gene_row = pan.index.get_loc(gene)
        
        # if new contig: resolve interval 
        if row.seqid != contig:
            pan.loc[interval, "posind"] = left_posind
            contig = row.seqid
            interval = []
            left_posind = None
            right_posind = None

        orthogroup = pan.iat[gene_row, orthogroup_col]

        if orthogroup in core: 

            # define position indicator at right side of interval
            # downstream means: interval is downstream of core orthogroup
            # --> core orthogroup is on minus strand
            downstream = row.strand == "-"
            right_posind = orthogroup + ("+" if downstream else "-")

            # assign right side position indicator to genes in interval
            # (left side position indicator would also be OK but not guaranteed
            # to be present)
            pan.loc[interval, "posind"] = right_posind

            # merge positions of left and right side position indicators
            if left_posind:
                left_position = positions.at[left_posind, "position"] 
                right_position = positions.at[right_posind, "position"] 
                if left_position != right_position:
                    has_right_position = positions.position == right_position
                    positions.loc[has_right_position, "position"] = \
                        left_position
                    
            if left_posind and intervals_file: 
                with open(intervals_file, "a") as handle:
                    genome = pan.iat[gene_row, genome_col]
                    if left_posind < right_posind:
                        handle.write(
                            ",".join([left_posind, right_posind, genome]))
                    else: 
                        handle.write(
                            ",".join([right_posind, left_posind, genome]))
                    handle.write("\n")

            interval = []
            left_posind = None
            right_posind = None

            # define position indicator at left side of new interval
            downstream = row.strand == "+"
            left_posind = orthogroup + ("+" if downstream else "-")

        else: 

            interval.append(gene)
            pan.iloc[gene_row, contig_col] = row.seqid
            pan.iloc[gene_row, strand_col] = row.strand
            pan.iloc[gene_row, start_col] = row.start
            pan.iloc[gene_row, end_col] = row.end

    return(positions, pan)
