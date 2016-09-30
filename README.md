#Output format


##chr    start	end	strand	outgroups	outgroup_chr positions  alleles strands score
	
- **chr** is the chromosome
- **start** 0-based start position on the chromosome for the site in the reference species genome
- **end** end position (end of range)
- **outgroups** are comma delimited list of outgroup species
- **outgroup_chr** are a comma delimited list of outroups species chromosomes for the aligned position (order as in the species listed in outgroups column)
- **positions** are the 0-based positions for the outgroups species 
- **alleles** are the nucleotides for the outgroup species alignment at that position in there genome
- **strands** are the strands for the nucleotides listed in alleles column
- **score** for the MAF alignment block

```
chr8    30876827        30876828        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8    30876827,27030444,4090372,6147243       T,T,T,T +,-,+,- 170727.0
```
