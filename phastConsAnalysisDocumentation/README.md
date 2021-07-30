# Summary of the Phast pipeline application

Included are the primary steps for applying the phastCons analysis as implemented for the N-fixing conserved element analysis. Primarilly RunPhastConsAnalysis.sh illustrates three primary steps of this analysis but would not trivially run verbatim in other working/computing environments (notably updating the indicated paths to the necessary software would be necessary).

We've utilized the v1_5 distribution of PHAST pipeline software (http://compgen.cshl.edu/phast/downloads.php) and a piece of dependency software, CLA Pack v 3.2, to compile the PHAST software binaries. It is specifically phyloFit and phastCons that we've used. See https://github.com/CshlSiepelLab/phast for more infromation. Both pieces of software are included as source code .tgz files in this directory (phast_v1_5.tgz and clapack.tgz) for completeness. 

Supporting tools from the maf tool set were also used as indicated in RunPhastConsAnalysis.s, available from https://github.com/dentearl/mafTools. Citation:
Genome Res. 2014 Dec;24(12):2077-89. doi: 10.1101/gr.174920.114. Epub 2014 Oct 1. Alignathon: a competitive assessment of whole-genome alignment methods. The information at the gitHub site should be fairly clear on how to get it working. 

# Step 1: Application of phyloFit to infer a general model of sequence evolution

We utilize a species tree and multiple slignmetn data set as input and obtain a model of sequence evolution with a usage as follows:

./phyloFit -t *tree-file* -o *prefix*.background.phylofit *maf-file*
  
This produces a model of sequence evolution with a file name <prefix>.background.phyloFit.mod, with the following format:
  
cat outside.background.phyloFit.mod 
ALPHABET: A C G T 
ORDER: 0
SUBST_MOD: REV
TRAINING_LNL: -127390897.232472
BACKGROUND: 0.291701 0.208880 0.208386 0.291033 
RATE_MAT:
  -0.948488    0.172353    0.478836    0.297300 
   0.240691   -1.068454    0.161162    0.666602 
   0.670279    0.161544   -1.073788    0.241965 
   0.297981    0.478432    0.173252   -0.949666 
TREE: (Vitis_vinifera_masked:0.0876077,((Eucalyptus_grandis_masked:0.239972,(Citrus_sinensis_masked:0.152403,(Theobroma_cacao_masked:0.117544,Arabidopsis_thaliana_masked:0.257447):0.0356378):0.023048):0.0207042,((Manihot_esculenta_masked:0.124674,(Populus_trichocarpa_masked:0.141391,(Linum_usitatissimum_masked:0.236434,Ricinus_communis_masked:0.105213):0.0217764):0.0123385):0.0567741,Medicago_truncatula_masked:0.230456):0.0195221):0.0876077);
  
Here the for loop over the variable G is for the species groups for which the analysis is being applied. 
  
Note: the tree input given has been the full species tree. This is because 1) the source species of Medicago is technically in the outside species set maf file, and 2) the species set in the tree will be correctly parsed down to the set consistently represented in the .maf data internally to phyloFit. This ensures the source species is always included from the multiple alignment data, even if not explicitly in every species set.
  
# Step 2: Estimate conserved and non-conserved sequence evolution models
  
The next step is to apply phastCons to the full set of multiple alignment data as follows:
  
./phastCons --estimate-trees *prefix.phastCons_v2.tree* *sorted-converted-.maf* *prefix.background.phylofit.mod*

Which produces two output model files, *prefix*.phastCons_v2.tree.cons.mod and *prefix*.phastCons_v2.tree.noncons.mod
  
Here the for loop over the variable G is again for each species group for which the analysis is being applied. Note also that in the .sh script lines for this step there are two commands applied to the initial .maf file for each species group analysis:
  
python ParseMAF.py *original-maf* *genome-.fai* > *converted maf*
mafSorter --maf *converted-maf* --seq Medicago_truncatula_masked.All > *converted-and-sorted-maf*
  
The ParseMAF.py script takes the original maf data and converts the Medicago source species region coordinates into coordinates relative to the full genome and with a single sequence name, "Medicago_truncatula_masked.All." This isnot by design, but was made necessary because the PhastCons software gave errors if multiple sequence namess appeared for the source species. Hence this was necessary to convert the multiple alignment data into a format that would include all of the data. In shortl this enabled training the phastCons conserved and non-conserved models on the fully genome-wide set of data, and not a subset of the data for one chromosome of the Medicago genome at a time. The ParseMAF.py script does this convesion with help from input .fai data, hence the inclusion of MtrunA17r5.0-20161119-ANR.fasta.fai here for utility.
  
The subsequent application of mafSorter is mainly to ensure that the regions are sorted in these new genome wide coordinates. 
  
This and the phyloFit step are the two most computationally intensive steps and depending on the number of species and volume of the multiple alignment data used as input can take times on the scale of days to weeks to complete.  
  
# Step 3: Calling conserved regions from the genome-wide models of conserved and non-conserved states
  
The two models learned from the genome-wide multiple alignment data (<prefix>.phastCons_v2.tree.cons.mod and <prefix>.phastCons_v2.tree.noncons.mod) are then finally used in a final step, which is to call conserved regions with phastCons using these inferred models. Per the mentioned concern about requiring only one sequence to be represented at a time, this is done for the multiply aligned regions on each Medicago v5 genome chromosome at a time. Hence the .sh will show a loop over both species groups (G) and chromosome (C). In this step regions are written out with the following usage with PhastCons:
  
 ./phastCons --most-conserved *path*/*group.chromosome*.bed -seq *chromosome* *maf-for-chromosome*  *conserved-model*,*non-conserved-model* > *output-.wig-file-of conservation-scores*
  
# Step 4: Selecting regions 5 bp or greater in size
  
 Finally, from the .bed files of conserved regions called for each Medicago genome chromsome, we select those regions 5 bp or larger and combine them into a single .bed file as shown in the .sh script.
