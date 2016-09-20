#Output format


##chr    start	end	strand	outgroups	outgroup_chr positions  alleles strands score
	
- **chr** is the chromosome
- **start** 0-based start position on the chromosome for the site in the reference species genome
- **end** end position 
- **outgroups** are comma delimited list of outgroup species
- **outgroup_chr** are a comma delimited list of outroups species chromosomes for the aligne position (order as in the species listed in outgroups column)
- **positions** are the 0-based positions for the outgroups species 
- **alleles** are the nucleotides for the outgroup species alignment at that position in there genome
- **strands** ae the strands for the nucleotides listed in alleles column
- **score** for the MAF alignment block


    chr8    30876827        30876828        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876827,27030444,4090372,6147243       T,T,T,T +,-,+,- 170727.0
    chr8    30876828        30876829        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876828,27030445,4090373,6147244       A,A,G,G +,-,+,- 170727.0
    chr8    30876829        30876830        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876829,27030446,4090374,6147245       T,T,T,T +,-,+,- 170727.0
    chr8    30876830        30876831        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876830,27030447,4090375,6147246       G,A,G,G +,-,+,- 170727.0
    chr8    30876831        30876832        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876831,27030448,4090376,6147247       G,C,G,G +,-,+,- 170727.0
    chr8    30876832        30876833        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876832,27030449,NA,6147248    T,C,-,A +,-,+,- 170727.0
    chr8    30876833        30876834        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876833,27030450,NA,6147249    A,A,-,A +,-,+,- 170727.0
    chr8    30876834        30876835        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876834,27030451,4090377,6147250       C,C,C,C +,-,+,- 170727.0
    chr8    30876835        30876836        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876835,27030452,4090378,6147251       A,A,A,A +,-,+,- 170727.0
    chr8    30876836        30876837        +       Greattit,Chicken,Zebrafinch,Flycatcher  chr8,chr8,chr8,chr8     30876836,27030453,4090379,6147252       T,T,C,T +,-,+,- 170727.0
