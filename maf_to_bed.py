from __future__ import print_function
import sys
import gzip
import argparse


def complement(seq):

    d = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', '-': '-', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}

    comp_seq = ''
    for n in seq:
        try:
            comp_seq += d[n]
        except KeyError:
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
    sites = {}
    positions = []
    strands = []

    gap_count = {spc: 0 for spc in spec}  # dictionary for holding the count of '-' characters in alignment block

    # print(len(aln_block[ref][5]))
    for pos in range(len(aln_block[ref][5])):

        # delayed printing
        if pos != 0 and aln_block[ref][5][pos] != '-':
            bed_line_str = [str(s) for s in bed_line]

            bed_line_str.append(','.join(species_lst))
            bed_line_str.append(','.join(chroms))
            bed_line_str.append(','.join(positions))
            bed_line_str.append(','.join([sites[si] for si in species_lst]))
            bed_line_str.append(','.join(strands))
            bed_line_str.append(score)

            print('\t'.join(bed_line_str))

            del bed_line_str[4:]
            del species_lst[:]
            del chroms[:]
            sites.clear()
            del positions[:]
            del strands[:]

            # print(bed_line)
            bed_line[1] += 1
            bed_line[2] += 1

        ins_rel_ref = False
        if aln_block[ref][5][pos] == '-':
                # gap_count[ref] += 1   # HB: used to keep ref coords correct after insertions relative to ref
                ins_rel_ref = True  # allows base extraction and pos update without chr  and species etc being updated

        for sp in spec:
            if sp not in species_lst:
                species_lst.append(sp)
                sites[sp] = ''

            # start.append(aln_block[sp][1] + site_num)
            if sp in aln_block.keys():
                sites[sp] += (aln_block[sp][5][pos])
                if ins_rel_ref is False:
                    chroms.append(aln_block[sp][0])
                    strands.append(aln_block[sp][3])

                if aln_block[sp][5][pos] == '-':
                    gap_count[sp] += 1
                    if ins_rel_ref is False:
                        positions.append('NA')
                else:
                    if ins_rel_ref is False:
                        positions.append(str(int(aln_block[sp][1]) + pos - gap_count[sp]))

            else:
                if ins_rel_ref is False:
                    sites.append('?')
                    chroms.append('?')
                    strands.append('?')
                    positions.append('?')

    # print final record in block
    bed_line_str = [str(s) for s in bed_line]

    bed_line_str.append(','.join(species_lst))
    bed_line_str.append(','.join(chroms))
    bed_line_str.append(','.join(positions))
    bed_line_str.append(','.join([sites[si] for si in species_lst]))
    bed_line_str.append(','.join(strands))
    bed_line_str.append(score)

    print('\t'.join(bed_line_str))

    del bed_line_str[4:]
    del species_lst[:]
    del chroms[:]
    sites.clear()
    del positions[:]
    del strands[:]

    # print(bed_line)
    bed_line[1] += 1
    bed_line[2] += 1


def main():

    parser = argparse.ArgumentParser(description="Convert a whole genome alignment in MAF format to BED format "
                                                 "following the coordinates of one of the species in the alignment ")
    parser.add_argument('-i', '--infile', dest='infile', required=True, help="Whole genome alignment in MAF format "
                                                                             "(Compressed)")
    parser.add_argument('-r', '--ref_species', dest='ref_species', required=True, help="Name of reference "
                                                                                       "species (as it appears "
                                                                                       "in the MAF file")
    parser.add_argument('-s', '--species_list', dest='species_list', required=True, help="Text file listing the "
                                                                                         "species names of the "
                                                                                         "species in the MAF file. This"
                                                                                         " file determines in which "
                                                                                         "order the species are in the "
                                                                                         "BED file")
    parser.add_argument('-c', '--chromosome', dest='ref_chrom', required=True, help="Specify which chromosome to"
                                                                                           "extract. This script only "
                                                                                           "extracts one to BED format "
                                                                                           "one chromosome in each run")

    args = parser.parse_args()

    ref_species = args.ref_species

    species_list = [i.rstrip() for i in open(args.species_list, 'r')]  # list of species to extract

    ref_chrom = args.ref_chrom

    with gzip.open(args.infile) as infile:
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

if __name__ == '__main__':
    main()


