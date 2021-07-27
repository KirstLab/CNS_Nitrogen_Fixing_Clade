[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)  
<!--[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXX) -->

<p align="center">
  <h2 align="center">Functional and comparative genomics reveals conserved noncoding sequences in the nitrogen-fixing clade</h2>
</p>

<!-- ABOUT THE PROJECT -->
## About this repository (in progress)

This repository contains the necessary code to replicate the results presented in the manuscript "Functional and comparative genomics reveals conserved noncoding sequences in the nitrogen-fixing" [[1]](#1).

## Prerequisites (in progress)

The code deposit in this repository integrates various software broadly applied in bioinformatics as well as in-home scripts. While the source code allows the reproduction of the results described in the paper, it is hard-coded for the genomes explored in the manuscript and tuned for execution in a [Slurm](https://slurm.schedmd.com/documentation.html) environment.

To generate the results presented in the publication, it is necessary the execution of multiple analytical steps. In this repository, the source code is organized in sessions the mimetic the mains analytical steps.

## Installation and execution (in progress)

### Download of genomes and whole-genome alignments

#### Step1 - Download of the genomes

The list of all genomes used in the manuscript is available on its supplemental file 1, including the web address to the download of each species. Importantly, you should download the unmasked version of the genomes.

#### Step2 - Genome masking using tantan

As mentioned, an unmasked version of the genomes is expected and the masking of repetitive elements is executed as part of the alignment pipeline, by applying "tantan" software [[2]](#2). The script *Whole-Genome_Alignments/tantan.sbatch* contains the code to mask the fasta files of all genomes.

#### Step3 - Pairwise alignments

The whole-genome alignments were implemented appling the [CNS pipeline](https://github.com/liangpingping/CNSpipeline). [[3]](#3)

## Usage (in progress)

## Reference
<a id="1">[1]</a> Pereira WJ, Knaack S, Conde D, Chakraborty S, Folk RA, Triozzi PM, Balmant KM, Dervinis C, Schmidt HW, Jean-Michel A, Roy S, Kirst M. [Functional and comparative genomics reveals conserved noncoding sequences in the nitrogen-fixing clade](), *BioRxiv*, 2021.

<a id="2">[2]</a> Frith, MC. [A new repeat-masking method enables specific detection of homologous sequences](https://academic.oup.com/nar/article/39/4/e23/1006710), *Nucleic Acids Research*, 2010.

<a id="3">[3]</a> Liang P. et al. [Single-Base Resolution Map of Evolutionary Constraints and Annotation of Conserved Elements across Major Grass Genomes]("https://academic.oup.com/gbe/article/10/2/473/4824837"), *Genome Biology and Evolution*, 2018.

<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3.0. See `LICENSE` for more information.
