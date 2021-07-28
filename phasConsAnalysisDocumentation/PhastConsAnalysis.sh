#Step 1: application of phyloFit. 

#point to executable and maf multiple alignment data directory. 
export exe=../../programs/PHAST_pipeline/phast/bin/phyloFit
export DATADIR=../../data/multiple_alignment_sequence_data/alignment_data/alignments_files/multiple_alignments

declare -A MAF
MAF[all_Nfix]=all_Nfix.sing.maf
MAF[outside]=outside.sing.maf

declare -A TREE
TREE[all_Nfix]=Nfix_subset.tree
TREE[outside]=outside_subset.tree

#output directory for fit models - variable for local analyses
OUTDIR1=../../results/sequence_analysis/phastCons/phyloFit_tree_models

rm commands_phyloFit.txt

for G in all_Nfix outside
do
    printf "${exe} -t all_genomes.tre -o ${OUTDIR1}/{G}.background.phyloFit ${DATADIR}/${MAF[$G]}\n" >> commands_PhyloFit.txt
done

#running multiple processes in parallel
cat commands_phyloFit.txt | xargs -d'\n' --max-procs=4 -I CMD bash -c CMD

#Step 2: application of phastCons to learn full model of conserved and diverged states using the -estimate tree option. 

if [ -e commandsPC.txt ]
then
    rm commandsPC.txt
fi

export exePC=../../programs/PHAST_pipeline/phast/bin/phastCons

FAI=MtrunA17r5.0-20161119-ANR.fasta.fai

for G in all_Nfix outside
do
	#convert input maf file in specific way: change the coodinates and naming of the separate Medicago chromosomes in order to represent the entire multiple alignment data with a single source sequence name, to prepare a model that is based on the genome wide set of multiple as much as possible. 
	#This is per the experience that otherwise it wasn't possible to obtain a fully genomewide model because multiple sequence names were producing an errors from the phastCons software.
	#The parseMAF.py script reorganizes the source species sequence names and coordiantes to represent everything in coodinates to the full Medicago genome, allowing the complete multiple alignment to be used.
	python ParseMAF.py ${datadir}/${MAF[$G]} ${FAI} > ${outdir}/training_maf/${G}.converted.maf 
	#with the converted sequence name and coordinate for Medicago implemented above in .coverted.maf, apply the mafTools sorter tool to make sure the blocks are coordinate-sorted. 
    mafSorter --maf ${outdir}/training_maf/${G}.converted.maf --seq Medicago_truncatula_masked.All > ${outdir}/training_maf/${G}.converted.sorted.maf
	#now list the commands to be issued to learn the "genomewide" models 
	printf "${exePC} --estimate-trees ${outdir}/models/${G}.phastCons_v2.tree ${outdir}/training_maf/${G}.converted.sorted.maf ${G}.background.phyloFit.mod > ${outdir}/model_training_wig_results/${G}.phastCons_v2.scores.wig\n" >> commandsPC.txt

done
#run the phastCons commands to learn genome-wide models with PhastCons
cat commandsPC.txt | xargs -d'\n' --max-procs=2 -I CMD bash -c CMD

#the output of each command will include an evolutionary model for a conserved and non-conserved state, both with prefixes ${outdir}/models/${G}.phastCons_v2.tree and file name endings of .cons.mod and .noncons.mod, respectively.

Step 3: Run PhastCons using the learned models to call conserved regions for each chromosome.

#define output directory for region calling from PhastCons
RDIR=../../results/sequence_analysis/phastCons/phastCons_called_regions

if [ -e commandsPC_perChrScoring.txt ]
then
	rm commandsPC_perChrScoring.txt
fi

for G in outside all_Nfix
do
	#loop over chromosomes
    for C in `cat ${FAI} | cut -f1 | sort -u`
    do
		#get .maf alignment blocks for this specific Medicago chromosome C using the maf tools extractor tool.
        ../../programs/maftools/mafTools/bin/mafExtractor --maf ${DATADIR}/${MAF[$G]} --seq Medicago_truncatula_masked.${C} --start 0 --stop 500000000 > ${OUTDIR}/per_chr_maf/${G}.${C}.sing.maf
		#make sure those blocks are coordinate sorted sorted 
        ../../programs/maftools/mafTools/bin/mafSorter --maf ${OUTDIR}/per_chr_maf/${G}.${C}.sing.maf --seq Medicago_truncatula_masked.${C} > ${OUTDIR}/per_chr_maf/${G}.${C}.sing.sorted.maf
		#List input files for region calling, including the extracted/sorted .maf file for multiple alignment blocks associated with a Medicago sequence on this chromosome C
        MAF=${OUTDIR}/per_chr_maf/${G}.${C}.sing.sorted.maf
        MA1=${OUTDIR}/models/${G}.phastCons_v2.tree.cons.mod
        MA2=${OUTDIR}/models/${G}.phastCons_v2.tree.noncons.mod
		#list appropriate PhastCons command
        printf "${exePC} --most-conserved ${RDIR}/${G}.${C}.bed --seqname ${C} ${MAF} ${MA1},${MA2} > ${OUTDIR}/genome_wide_model_results/${G}.${C}.scores.wig\n" >> commandsPC_perChrScoring.txt
    done
done

#run full set of phastCons commands, will complete fairly quickly because models are given as input. 
cat commandsPC_perChrScoring.txt | xargs -d'\n' --max-procs=7 -I CMD bash -c CMD

Step 4 - Filter regions for size. 

#for each set of regions called per chromosome, organize a summary list for all chromosomes selecting those that are 5 bp in size.
for G in all_Nfix outside
do
    cat ${RDIR}/${G}.Mtrun*.bed | awk '$3-$2>=5{print}' > ${DATADIR}/${G}.all.bed
done

