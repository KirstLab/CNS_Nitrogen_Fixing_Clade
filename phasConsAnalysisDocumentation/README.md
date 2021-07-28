# Summary of the Phast pipeline application

Included are the primary steps for applying the phastCons analysis as implemented for the N-fixing conserved element analysis. Primarilly RunPhastConsAnalysis.sh illustrates three primary steps of this analysis. 

Here we've utilized the v1_5 distribution of PHAST pipeline software (http://compgen.cshl.edu/phast/downloads.php) and a piece of dependency software CLA Pack v 3.2 to compile the PHAST software binaries, specifically phyloFit and phastCons that  we've utilized. Also see https://github.com/CshlSiepelLab/phast for more infromation. Both are included as .tgz files in this directory for completeness. 

Supporting tools from the maf tool set were also used as indicated in RunPhastConsAnalysis.s, available from https://github.com/dentearl/mafTools. Citation:
Genome Res. 2014 Dec;24(12):2077-89. doi: 10.1101/gr.174920.114. Epub 2014 Oct 1. Alignathon: a competitive assessment of whole-genome alignment methods.

# Step 1: Application of phyloFit to infer a general model of sequence evolution

