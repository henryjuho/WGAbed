# WGAbed - a package for handling whole genome alignments

This packages handles the conversion of a whole genome alignment file in MAF format to a reference ordered BED format containing all the information in the MAF file but in a more accessible format. This respository contains scripts for creating a 'whole genome alignment BED' or 'WGAbed' file as well as scripts for downstream manipulation of the file.

## Creating the whole genome alignment bed file

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

## polarise_wga_ref_indels.py

This script returns off wga bed lines that contain the specified INDEL type.

### Usage

```
$ ./polarise_wga_ref_indels.py -h
usage: polarise_wga_ref_indels.py [-h] -indel_type
                                  {non_indel,ambig_indel,deletion,insertion}

optional arguments:
  -h, --help            show this help message and exit
  -indel_type {non_indel,ambig_indel,deletion,insertion}
                        type of indel to output
```

### Example and output

To extract all ref specific insertions:

```
$ ./maf_to_bed.py -i data/test.maf.gz -r Zebrafinch -c chr8 | sort -k1,1 -k2,2n | ./polarise_wga_ref_indels.py -indel_type insertion
chr8    4089946 4089971 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4089946,27030055,6146813,30876430       TTTGTATGAATTCTTTTGAATACCT,T------------------------,T------------------------,T------------------------   +,-,-,+ 170727.0
```

To extract all ref specific deletions:

```
$ ./maf_to_bed.py -i data/test.maf.gz -r Zebrafinch -c chr8 | sort -k1,1 -k2,2n | ./polarise_wga_ref_indels.py -indel_type deletion
chr8    4090041 4090042 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090041,27030121,6146884,30876501       T-,TG,TG,TG     +,-,-,+ 170727.0
chr8    4090219 4090220 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090219,27030264,6147064,30876680       G--,GAG,GAG,GAG +,-,-,+ 170727.0
chr8    4090314 4090315 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090314,27030361,6147161,30876770       C-,TC,TT,TT     +,-,-,+ 170727.0
chr8    4090374 4090375 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090374,27030423,6147222,30876831       G--,CCA,GAA,GTA +,-,-,+ 170727.0
chr8    4090379 4090380 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090379,27030430,6147229,30876838       G-------,GTGCTAAT,GTGTTAAT,GTGTTAAT       +,-,-,+ 170727.0
```

To extract all ambiguous INDELs:

```
$ ./maf_to_bed.py -i data/test.maf.gz -r Zebrafinch -c chr8 | sort -k1,1 -k2,2n | ./polarise_wga_ref_indels.py -indel_type ambig_indel
chr8    4089892 4089896 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4089892,27030005,6146759,30876376       TAGT,A---,CAGT,TAGT     +,-,-,+ 170727.0
chr8    4089943 4089945 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4089943,27030053,6146810,30876427       AC,A-,AC,AC     +,-,-,+ 170727.0
chr8    4090000 4090003 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090000,27030085,6146843,30876460       TTT,C--,ttt,TTT +,-,-,+ 170727.0
chr8    4090022 4090026 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090022,27030105,6146865,30876482       ATTG,A---,TTTG,ATTG     +,-,-,+ 170727.0
chr8    4090076 4090077 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090076,27030157,6146920,30876537       A-,C-,AA,A-     +,-,-,+ 170727.0
chr8    4090080 4090100 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090080,27030161,6146925,30876541       TAGTATTTCTGGATTAATTT,T-------------------,TGGtatttctggatttattt,TTCTATTTCTGGATTAATTT       +,-,-,+ 170727.0
chr8    4090120 4090125 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090120,27030182,6146965,30876581       TGTGT,T----,tgtga,CATGA +,-,-,+ 170727.0
chr8    4090187 4090201 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090187,27030245,6147032,30876648       CTGCTCTGAAAGTT,C-------------,CTGCTCTGAAAGTC,CTGCTCTAAAAGTC       +,-,-,+ 170727.0
chr8    4090269 4090277 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090269,27030316,6147116,30876732       AAAGAAGA,AAAGAAGG,aaagaaga,A-------       +,-,-,+ 170727.0
chr8    4090361 4090362 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090361,27030409,6147209,30876818       G-,AT,G-,C-     +,-,-,+ 170727.0
chr8    4090407 4090413 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090407,27030465,6147264,30876873       GTGTCA,G-----,GCTTGG,GTATCA     +,-,-,+   170727.0
chr8    4090424 4090425 +       Zebrafinch,Chicken,Flycatcher,Greattit  chr8,chr8,chr8,chr8     4090424,27030477,6147281,30876890       A-,GA,A-,AA     +,-,-,+ 170727.0
```

