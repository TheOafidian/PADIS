# PADIS: pangenome-based discovery of insertion sequences

PADIS is a tool for the discovery of insertion sequences in prokaryotic genomes in a database-independent manner. In a nutshell, PADIS identifies orthogroups (gene families) present in more than one location within or across genomes and then checks them for the presence of terminal inverted repeats (TIRs). The tool requires a set of genome assemblies (.fna files), gene annotation for these assemblies (.gff files) and a pangenome in [SCARAP](https://github.com/SWittouck/SCARAP) format.

PADIS is still very much in active development. This means that the interface and output may change rapidly and drastically, that there may be bugs and that some types of input may not (yet) work. Even the name may still change. Feel free to post [issues](https://github.com/swittouck/PADIS/issues).

## Installation

When PADIS is more mature, an official release will be uploaded to PyPI and a recipe to bioconda.

For now, you can install PADIS directly from GitHub with pip. There are only Python dependencies.

    pip install git+https://github.com/swittouck/PADIS

## How to run

Given a set of assemblies (.fna files) in an `assemblies` folder, first perform gene prediction, e.g. with Prodigal:

    mkdir genes
    for dir in gffs faas logs; do mkdir genes/$dir; done
    for assembly_file in assemblies/*.fna; do
      genome=$(basename $assembly_file .fna)
      prodigal \
        -i $assembly_file \
        -f gff \
        -o genes/gffs/${genome}.gff \
        -a genes/faas/${genome}.faa \
        2> genes/logs/${genome}.txt
    done

Then infer the pangenome, e.g. with SCARAP, using 8 threads:

    scarap pan genes/faas pangenome -t 8

Finally, run PADIS as follows:

    padis assemblies genes/gffs pangenome/pangenome.tsv padis -t 8

This produces the following files in the `padis` folder:

* `accessory_genes.csv`: gene table with a bunch of columns; the most important one is "position", which contains a position identifier for the gene
* `accessory_orthogroups.csv`: most important output: all information relevant to the prediction of the orthogroups as IS orthogroups (i.e., transposases or regulatory genes):
    * orthogroup: identifier of the orthogroup
    * genes: number of genes in the orthogroup
    * genomes: number of genomes the orthogroup is present in
    * located: number of genes where a position could be determined
    * positions: number of unique positions the orthogroup is present in across the genomes
    * status: either "potential IS" or a reason why the orthogroup was not identified as a potential IS
    * length: length of the homologous region around a representative gene of the orthogroup (determined by aligning against another gene of the orthogroup at a different location, including flanking regions)
    * tir_score: alignment score of the 30 first nucleotides of the homologous region against the reverse complement of the last 30 nucleotides
    * tir_nscore: normalized version of the tir_score using 100 alignments of shuffled versions of the terminal sequences: `nscore = (score - mean(random_scores)) / sd(random_scores)`
    * tir_length: length of the TIR
* `summary.csv`: orthogroup counts for the different "status" options

## Licence

PADIS is free software, licensed under [GPLv3](https://github.com/SWittouck/padis/blob/master/LICENSE).
