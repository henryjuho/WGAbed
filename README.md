# Conversion of whole genome alignment file in MAF format to a reference ordered BED format 

##Example usage

    $ ./maf_to_bed.py -i data/test.maf.gz -r Greattit -s data/species.txt -c chr8 | sort -k1,1 -k2,2n | gzip -c > data/test.wga.bed.gz



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


##Example aligned site in the BED format
```
chr8    30876827        30876828        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8    30876827,27030444,4090372,6147243       T,T,T,T +,-,+,- 170727.0
```

#Utility scripts
##wga_bed_indels.py

This script takes the piped output of ```maf_to_bed.py``` and extracts the INDELs. It is possible to specify a number of filtering options, see usage below.
 
###Usage

```
$ ./wga_bed_indels.py -h
usage: wga_bed_indels.py [-h] [-max_length MAX_LENGTH]
                         [-min_coverage MIN_COVERAGE] [-ref_specific]

optional arguments:
  -h, --help            show this help message and exit
  -max_length MAX_LENGTH
                        Maximum INDEL length to extract
  -min_coverage MIN_COVERAGE
                        Minimum species coverage
  -ref_specific         Restricts output to INDELs only found in reference
                        species, and conserved in non reference species
```

###Example and output

```
$ ./maf_to_bed.py -i data/test.maf.gz -r Zebrafinch -c chr8 | sort -k1,1 -k2,2n | ./wga_bed_indels.py -min_coverage 4 -max_length 10 -ref_specific 
```

```
chr8	4090041	4090042	+	Zebrafinch,Chicken,Flycatcher,Greattit	chr8,chr8,chr8,chr8	4090041,27030121,6146884,30876501	T-,TG,TG,TG	+,-,-,+	170727.0
chr8	4090219	4090220	+	Zebrafinch,Chicken,Flycatcher,Greattit	chr8,chr8,chr8,chr8	4090219,27030264,6147064,30876680	G--,GAG,GAG,GAG	+,-,-,+	170727.0
chr8	4090314	4090315	+	Zebrafinch,Chicken,Flycatcher,Greattit	chr8,chr8,chr8,chr8	4090314,27030361,6147161,30876770	C-,TC,TT,TT	+,-,-,+	170727.0
chr8	4090374	4090375	+	Zebrafinch,Chicken,Flycatcher,Greattit	chr8,chr8,chr8,chr8	4090374,27030423,6147222,30876831	G--,CCA,GAA,GTA	+,-,-,+	170727.0
chr8	4090379	4090380	+	Zebrafinch,Chicken,Flycatcher,Greattit	chr8,chr8,chr8,chr8	4090379,27030430,6147229,30876838	G-------,GTGCTAAT,GTGTTAAT,GTGTTAAT	+,-,-,+	170727.0
```