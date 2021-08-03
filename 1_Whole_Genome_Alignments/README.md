## Download of genomes and masking of repetitive elements

The list of all genomes used in the manuscript is available on its Supplemental File 1, including the web address to the download of each species. Importantly, the unmasked version of the genomes are required. When using LAST to promote the alignments, it is recommended that the masking of repetitive elements is executed by applying the tantan software [[1]](#1). In the file *tantan.sbatch* a piece of code to remove softmasking is available as well as the code to execute tantan.

### Usage
```shell
sbatch tantan.sbatch
```

## Pairwise alignments

The pairwise whole-genome alignments were implemented using a modified version of the [CNS pipeline](https://github.com/liangpingping/CNSpipeline) [[2]](#2). The changes were made aiming to allow the execution using multiple arrays and the modified scripts are available in the folder CNSpipeline_modified/. The scripts below assume that the CNSpipeline_modified/ files are accessible in your path.

### Usage

First, it is necessary to prepare the reference and target genomes

```shell
# Prepare the reference genome
sbatch whole_genome_align_prep_reference.sbatch
# Prepare the target genomes
sbatch whole_genome_align_prep_other_genomes.sbatch
```

Then, execute the pairwise whole-genome alignment of each genome against the reference. The steps below must be executed for each genome. Users need to replace <genome_species> by the name of the genome to be aligned.

```shell
# Executes the pairwise whole-genome alignment
sbatch genome_alig_part1.sbatch
# After part 1 is concluded
sbatch genome_alig_part2.sbatch
```

After the execution of the pairwise alignments of all genomes.
```shell
sbatch genome_alignments_chaining.sbatch
```

## Generating multiple alignments using ROAST

The pairwise-alignment generated are stored in the reference/tba/ folder and the next steps assume execution in this directory (or you can copy the maf files to your directory of choice).

[ROAST](https://www.bx.psu.edu/~cathy/toast-roast.tmp/README.toast-roast.html) joins the pairwise alignments using a phylogenetic tree as guide. ROAST version used in this study is available at: https://www.bx.psu.edu/miller_lab/dist/multiz-tba.012109.tar.gz, last accessed in August 2021.

We found that the execution works better when ROAST is used to export its commands lines as a shell script, that is executed in opposition to executing ROAST directly.





## References
<a id="1">[1]</a> Frith, MC. [A new repeat-masking method enables specific detection of homologous sequences](https://academic.oup.com/nar/article/39/4/e23/1006710), *Nucleic Acids Research*, 2010.

<a id="2">[2]</a> Liang P. et al. [Single-Base Resolution Map of Evolutionary Constraints and Annotation of Conserved Elements across Major Grass Genomes]("https://academic.oup.com/gbe/article/10/2/473/4824837"), *Genome Biology and Evolution*, 2018.
