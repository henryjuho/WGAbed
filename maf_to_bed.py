from __future__ import print_function
import sys
import gzip


def complement(seq):

    comp_seq = ''
    for n in seq.upper():
        if n == 'A':
            comp_seq += 'T'
        elif n == 'T':
            comp_seq += 'A'
        elif n == 'G':
            comp_seq += 'C'
        elif n == 'C':
            comp_seq += 'G'
        elif n == '-':
            comp_seq += '-'
        elif n == 'N':
            comp_seq += 'N'
        else:
            sys.exit('Unrecognised aligment character %s' % n)

    return comp_seq


def revcomp(seq):

    rev_seq = seq[::-1]
    rc_seq = complement(rev_seq)

    return rc_seq


def revpos(pos, seq_len, block_len):

    reverse_pos = int(seq_len) - int(pos) - block_len

    return reverse_pos


def create_bed_records(aln_block, spec, ref, score):

    bed_line = [aln_block[ref][0], int(aln_block[ref][1]), int(aln_block[ref][1]) + 1, aln_block[ref][3]]

    species_lst = []
    chroms = []
    sites = []
    positions = []
    strands = []

    gap_count = {spc: 0 for spc in spec} # dictionary for holding the count of '-' characters in alignment block

    # print(len(aln_block[ref][5]))
    for pos in range(len(aln_block[ref][5])):
        if aln_block[ref][5][pos] == '-':
            gap_count[ref] += 1
            continue
        for sp in spec:
            species_lst.append(sp)
            # start.append(aln_block[sp][1] + site_num)
            if sp in aln_block.keys():
                sites.append(aln_block[sp][5][pos])
                chroms.append(aln_block[sp][0])
                strands.append(aln_block[sp][3])

                if aln_block[sp][5][pos] == '-':
                    positions.append('NA')
                    gap_count[sp] += 1
                else:
                    positions.append(str(int(aln_block[sp][1]) + pos - gap_count[sp]))

            else:
                sites.append('?')
                chroms.append('?')
                strands.append('?')

        bed_line_str = [str(s) for s in bed_line]

        bed_line_str.append(','.join(species_lst))
        bed_line_str.append(','.join(chroms))
        bed_line_str.append(','.join(positions))
        bed_line_str.append(','.join(sites))
        bed_line_str.append(','.join(strands))
        bed_line_str.append(score)

        print('\t'.join(bed_line_str))

        del bed_line_str[4:]
        del species_lst[:]
        del chroms[:]
        del sites[:]
        del positions[:]
        del strands[:]

        # print(bed_line)
        bed_line[1] += 1
        bed_line[2] += 1


ref_species = sys.argv[2]

species_list = [i.rstrip() for i in open(sys.argv[3], 'r')]  # list of species to extract

ref_chrom = sys.argv[4]

with gzip.open(sys.argv[1]) as infile:
    # align_block = {}
    for line in infile:
        if line.startswith('#'):
            continue
        elif line.startswith('a'):
            align_score = line.split()[1].split('=')[1]
            align_block = {}
            continue
        elif line.startswith('s'):
            species = line.split()
            align_block[species[1].split('.')[0]] = [species[1].split('.')[1]] + species[2:]

        else:
            # block_start = 0
            if ref_species not in align_block.keys():
                align_block = {}
                continue
            else:
                ref_seq_entry = align_block[ref_species]
                chrom = ref_seq_entry[0]
                if chrom != ref_chrom:
                    continue
                ref_strand = ref_seq_entry[3]
                # print(align_block)
                if ref_strand == '-':
                    for s in align_block.keys():
                        align_block[s][1] = revpos(align_block[s][1], align_block[s][4], int(align_block[s][2]))
                        align_block[s][5] = revcomp(align_block[s][5])
                        if align_block[s][3] == '-':
                            align_block[s][3] = '+'
                        else:
                            align_block[s][3] = '-'

                create_bed_records(align_block, species_list, ref_species, align_score)



##TODO check if the reference genomome is covered by more than one alignment in the MAF file


