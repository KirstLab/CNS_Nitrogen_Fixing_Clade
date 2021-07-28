import sys
import math
import os

f = open(sys.argv[1],'r')

fai = open(sys.argv[2],'r')

chrDict={}
#initialize 
count=0
for line in fai:
	info=line.strip().split('\t');
	#get effective starting coordinate of this chromosome in the full genome-wide seuqence coordinate space, i.e, that spanning all chromosomes
	chrDict[info[0]]=count+1
	#increment count to be the starting coordinate for the next-included chromosome
	count=count+int(info[1])


#for c in chrDict:
#	print c
#	print chrDict[c]

for line in f:
	parsed=line.strip().split();
	#check if line is associated with a Medicago sequence as a source species 
	if parsed and parsed[0] == str('s') and parsed[1].strip().split('.')[0] == str('Medicago_truncatula_masked'):
		#if so get chromosome name of original sequence from Medicago 
		chromosome=parsed[1].strip().split('.')[1]
		#write out same information in new "genome-wide" coordinates based on the fai information for the Medicago genome. Give this sequence a new sequence ".All" sequence name, and coordinates relative to the whole genome. 
		print parsed[0]+' Medicago_truncatula_masked.All '+str(int(parsed[2])+int(chrDict[chromosome])-1)+' '+parsed[3]+' '+parsed[4]+' 436900332 '+parsed[6]
		#donen becaues the "sequnece name" has to be uniform in the phastCons analysis, and we want a model that is as genomewide as possible
	else:
		#otherwise write out information unchanged. 
		sys.stdout.write(line)
