# Conversion of whole genome alignment file in MAF format to a reference ordered BED format 

##Example usage

    $ python maf_to_bed.py -i data/test.maf.gz -r Greattit -s data/species.txt -c chr8 | sort -k1,1 -k2,2n | gzip -c > data/test.wga.bed.gz



#Output format


##chr    start	end	strand	outgroups	outgroup_chr positions  alleles strands score
	
- **chr** is the chromosome in the reference species genome
- **start** 0-based start position on the chromosome for the site in the reference species genome
- **end** end position (end of range) in the reference species genome
- **outgroups** are comma delimited list of outgroup species (Order determined based on species list input file)
- **outgroup_chr** are a comma delimited list of outroups species chromosomes for the aligned position 
(order as in the species listed in outgroups column)
- **positions** are the 0-based positions for the outgroups species. A '?' is used when a non-reference species is not 
present within a MAF block and '-' for INDELS in the alignment.
- **alleles** are the nucleotides for the outgroup species alignment at that position in there genome. A '?' is used when a non-reference species is not present within a MAF block
- **strands** are the strands for the nucleotides listed in alleles column
- **score** for the MAF alignment block


##Example aligne site in the BED format
```
chr8    30876827        30876828        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8    30876827,27030444,4090372,6147243       T,T,T,T +,-,+,- 170727.0
```
