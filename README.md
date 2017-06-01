# WGAbed - a package for handling whole genome alignments

This packages handles the conversion of a whole genome alignment file in MAF format to a reference ordered BED format containing all the information in the MAF file but in a more accessible format. This respository contains scripts for creating a 'whole genome alignment BED' or 'WGAbed' file as well as scripts for downstream manipulation of the file.

## Creating the whole genome alignement bed file

A WGAbed file can be created from a MAF file with the script ```maf_to_bed.py```. Script must be passed the species which you would like the BED coordinates to follow with the ```-r``` flag. Additionally the chromosome must be specified with ```-c```. The output can then be piped to sort and bgzip to produced a sorted, indexable WGAbed file:

```
$ ./maf_to_bed.py -i data/test.maf.gz -r Greattit -s data/species.txt -c chr8 | sort -k1,1 -k2,2n | bgzip -c > data/test.wga.bed.gz
```

This file can then be indexed with tabix:

```
$ tabix -pbed data/test.wga.bed.gz
```

## WGAbed format


## chr    start	end	strand	outgroups	outgroup_chr positions  alleles strands score
	
- **chr** is the chromosome in the reference species genome
- **start** 0-based start position on the chromosome for the site in the reference species genome
- **strand** the strand for the reference species sequence
- **end** end position (end of range) in the reference species genome
- **outgroups** are comma delimited list of outgroup species (Order determined based on species list input file)
- **outgroup_chr** are a comma delimited list of outroups species chromosomes for the aligned position 
(order as in the species listed in outgroups column)
- **positions** are the 0-based positions for the outgroups species. A '?' is used when a non-reference species is not 
present within a MAF block and '-' for INDELS in the alignment.
- **alleles** are the nucleotides for the outgroup species alignment at that position in their genome. A '?' is used when a non-reference species is not present within a MAF block
- **strands** are the strands for the nucleotides listed in alleles column
- **score** for the MAF alignment block


## Example aligned site in the BED format
```
chr8    30876827        30876828        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8    30876827,27030444,4090372,6147243       T,T,T,T +,-,+,- 170727.0
```

# Utility scripts
## wga_bed_indels.py

This script takes the piped output of ```maf_to_bed.py``` and extracts the INDELs. It is possible to specify a number of filtering options, see usage below.
 
### Usage

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

### Example and output

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

## non_ref_intersect.py

This script takes a piped wga.bed file and extracts a subset of sets specified by an associated bed file with coordinates from a non-reference species in the alignment. NOTE output is still written as a reference ordered bed.
 
### Usage

```
$ ./non_ref_intersect.py -h
usage: non_ref_intersect.py [-h] -b BED_NON_REF -q QUERY_SPECIES -c CHROMOSOME

optional arguments:
  -h, --help            show this help message and exit
  -b BED_NON_REF, --bed_non_ref BED_NON_REF
                        .bed file with non refernce spp coordinates to
                        intersect with wga.bed file
  -q QUERY_SPECIES, --query_species QUERY_SPECIES
                        species that the query bed file has coordinates for
  -c CHROMOSOME, --chromosome CHROMOSOME
                        target chromosome

```

### Example and output

```
$ cat data/test_indel_data.bed | ./non_ref_intersect.py -b data/test_query.bed.gz -q Flycatcher -c chr8
```

```
chr8    30876342        30876343        +       Greattit,Chicken,Flycatcher,Zebrafinch  chr8,chr8,chr8,chr8     30876342,27029971,6146725,4089858       A,a,A,A +,-,-,+ 170727.0
chr8    30876343        30876344        +       Greattit,Chicken,Flycatcher,Zebrafinch  chr8,chr8,chr8,chr8     30876343,27029972,6146726,4089859       A,a,A,A +,-,-,+ 170727.0
chr8    30876344        30876345        +       Greattit,Chicken,Flycatcher,Zebrafinch  chr8,chr8,chr8,chr8     30876344,27029973,6146727,4089860       A,t,A,A +,-,-,+ 170727.0
chr8    30876345        30876346        +       Greattit,Chicken,Flycatcher,Zebrafinch  chr8,chr8,chr8,chr8     30876345,27029974,6146728,4089861       A,a,A,A +,-,-,+ 170727.0
chr8    30876346        30876347        +       Greattit,Chicken,Flycatcher,Zebrafinch  chr8,chr8,chr8,chr8     30876346,27029975,6146729,4089862       T,t,T,T +,-,-,+ 170727.0

```
## maf_extract_ref_chr.py

This script extracts all blocks containing a given reference chromosome.

### Usage

```
$ ./maf_extract_ref_chr.py -h
usage: maf_extract_ref_chr.py [-h] -c CHROMOSOME [-H]

Extract a reference chromosome from a WGA MAF file

optional arguments:
  -h, --help            show this help message and exit
  -c CHROMOSOME, --chromosome CHROMOSOME
                        Specify which chromosome to extract
  -H, --header          If specified will print header in output
```

### Example and output

```
$ zgrep -v ^6 data/test.maf.gz | ./maf_extract_ref_chr.py -c chr8
a score=170727.0
s Chicken.chr8     2932491 551 + 29963013 CAATGGAATTTAGTAGATAAGACTCAGAATCAGATTCTGGGTTTCAAAGTTGAATG-----CAAGGAAGGAAAAAATCCCCACAAGTATATTAGCACAATGTGGTATAAATCCATCATCTATAATGTACTGGATTGCTTTTAACTCTCAGTCTGAAACCAAATCGATTTCAACAGTTGGTGCTGGACAGAGAACAATTACTAACCTTCTTTTTTCTCTTTCCTGCTGCTTAACAGGGCTACAAATAACAACAGAAAATTACTCTCAAAGATCTACCCATGT-------------GAAACCAAAGTAACTTAAATGTTCAAATATTATTTGGATTGCAATGTGTAATGTCTCTCAAAT----AGACAAGAAAAATATGGGCAA-------------------AATA-GTCTATAGTCAGTTTAAGAGTCCCTTCTACATTTGCAATCATGTAAAATGAT---TTAGGTATCAATTTCAAAAC--GTTCAAATACAGGAAAGAAAGTTCCTCAGT------------------------AT-TATTAACATATTGTAGAACATTCGCCATCATCCTGGCATTTCAAAACG---TTAAGTTTAATGAGCTAGCCAAGTAGttatatatt
s Zebrafinch.chr8 23902959 610 - 27993427 CACTGGGAGATGCTGAAAGAGCCCCGGGATGGGAGTGTGGGGT-TAAAGCTGAGTGTGACACAGGGCAGGCCACGAGCCCCTGAGTCCC-------CAGTG--CCACAGATCCAGC-CCTGCAATGCCCTGCCTGCCTTTCACCTCCCAGCCTGCAGCCAAACC-GTGTCAACACTTTGTGCTGGACACAGACCAATTACTAATCTTCTTTCTCTTCTTTCCTGCTGCTTAACAGGCCTACAAATAATGCTTTAAAATTA--CTCTGAGACCTACCCACATAACTTTCAGAGCAGAAACCAAAGTAACTTAAATTATCAAAAATCATTTCTATTTCAAGGTGAAATGTCTTTGAAATACACAGAAAATAAAAATCCTGGTTAAAATTAATCCAGAAATACTAATG-TATTGAGGTCACTGCAAGGATCCTATCTGGATTTG-AATCACGTAAAAGAATCAATTAAGAGTCAAATTTAAAAGAAATTCAAAAACAGGAAAGATAATTTCTCAGTAGGTATTCAAAAGAATTCATACAAATGTATTTAGGTGTTGCAGAACACTTACCATCACTCTAGCACTACCAAACAACTATAAATCTATTTAACTAACCTGTTAAGTGGATTTT
s Flycatcher.chr8 25953491 600 + 32100816 CACTGGAAAACGCTAAAAAAGGCTCAGGATGGGGCCATGGGGG-TGGAGCTGAGTGCCAAGCAGGGCAGGTGAGGAGCCCCTGGCGTCCATTAACACAATGTTCCACAGATCCATC-CCTGTAATGCACTGCCTGCCTTTTACCTCGCAGCCTGAAGCCAAATCAATTTCAACAGTTTGTGCTGGACACAGACCAATTACTaatcttctttctcctctttcctgctgcttAACAGGGCTACAAATAATGGTTTAAAATTACTCTCTGAGACCTACCCACATGACTTTCAGAGCAGAAACCAAAGTAACTTGAATGATCAAAAATCGTTTATACTTCAAGGTGAAATGTCtttgaaattcacagaaaagaaaaatcctggctaaaataaatccagaaataCCAACGTTTTTGTGGTGTCTGCAAGGATCCTATCTGGATTTGCAATCATGTGAAAGAATCAAATAGGagtcaaatttaaaataaattcaaaaacaGGAAAGGGAACTTCTCAGT------------------------ATGTATTCAGGTGTTGCAGAACACTTACCATCATTCTATCACTGCCAAACAACTGTAAATCTATTTAACCAAACAATTGAGTGGATTTT
s Greattit.chr8     447231 593 - 31324166 CACTGGAATATGCTAAAAAAGAATCAGGATGGGACCATGGGGGTTAAAGGTGAGTGTGATACAGGGAAGATGAAGAGCCCCTGGTGTCGATTAACACAATGTACCATAGATCCATC-GCTCTAATGCATTGTCTGGCTTTTACCTCGCAGCCTGAAGCCAAATCAATTTCAACAGTTTGTGCTGGACAGACACCAATTACAAA-------TTACCTCTTTCCTGCTGCTTAACAGGGCTACAAAAAATGGTTTAAAATTACTCTCTAAGACCTACCCACATGACTTTTAGAGCAGAAACCAACGTAACTTAAATGATCAAAAACTGTTTCTATTTCAAGGTGAAATGTCTTGGAAATTCATGGAAAAGAAAAATCCATGCTAAAATTAATCCAGAAATAGAAATG-TTTTCTGGTCACTGTAAGGAACTTAACTGGATTTGCAATCATAAAAAAGAATCAATTAGGAGTCAAATTTAAAATAAATTTAAAAATAGGAAAAGGAACTTCTCAGT------------------------ATGTATTCAGGTATTGCAGAACACTTAACATCATTCTGGCACTGCCAAACAACTATACGTCTATTTAACTAACTAACTAAGTGGATTTT
```
