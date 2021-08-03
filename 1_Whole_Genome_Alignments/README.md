## Installation and execution (in progress)

### Download of genomes and whole-genome alignments

#### Step1 - Download of the genomes

The list of all genomes used in the manuscript is available on its supplemental file 1, including the web address to the download of each species. Importantly, you should download the unmasked version of the genomes.

#### Step2 - Genome masking using tantan

As mentioned, an unmasked version of the genomes is expected and the masking of repetitive elements is executed as part of the alignment pipeline, by applying "tantan" software [[2]](#2). The script *Whole-Genome_Alignments/tantan.sbatch* contains the code to mask the fasta files of all genomes.

#### Step3 - Pairwise alignments

The whole-genome alignments were implemented appling the [CNS pipeline](https://github.com/liangpingping/CNSpipeline). [[3]](#3)

## Usage (in progress)


## References
<a id="2">[2]</a> Frith, MC. [A new repeat-masking method enables specific detection of homologous sequences](https://academic.oup.com/nar/article/39/4/e23/1006710), *Nucleic Acids Research*, 2010.

<a id="3">[3]</a> Liang P. et al. [Single-Base Resolution Map of Evolutionary Constraints and Annotation of Conserved Elements across Major Grass Genomes]("https://academic.oup.com/gbe/article/10/2/473/4824837"), *Genome Biology and Evolution*, 2018.
