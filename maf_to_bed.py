#!/usr/bin/env python

from __future__ import print_function
import sys
import gzip
import argparse
import subprocess


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
    ref_start = int(aln_block[ref][1])

    species_lst = []
    chroms = []
    sites = {}
    positions = []
    strands = []

    gap_count = {spc: 0 for spc in spec}  # dictionary for holding the count of '-' characters in alignment block
    start_print = False
    # print(len(aln_block[ref][5]))
    for pos in range(len(aln_block[ref][5])):

        # delayed printing
        if start_print is True and '-' not in [aln_block[x][5][pos] for x in aln_block.keys()]:
            ref_end = ref_start + (len(sites[ref]) - sites[ref].count('-'))
            bed_line_str = [str(s) for s in [bed_line[0], ref_start, ref_end, bed_line[3]]]

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

        indel = False
        # catches indels and allows bases to be appended to previous site instead of constructing new bed line
        if '-' in [aln_block[x][5][pos] for x in aln_block.keys()]:
            indel = True

        # suppress print for indels that occur at start of block, due to inability to assign genomic start positions
        if indel is False:
            start_print = True
        if start_print is False:
            # updates gap count prior to skipping site
            for sp in spec:
                if sp in aln_block.keys():
                    if aln_block[sp][5][pos] == '-':
                        gap_count[sp] += 1
            continue

        for sp in spec:
            if sp not in species_lst:
                species_lst.append(sp)
                sites[sp] = ''

            # start.append(aln_block[sp][1] + site_num)
            if sp in aln_block.keys():
                sites[sp] += (aln_block[sp][5][pos])
                if indel is False:
                    chroms.append(aln_block[sp][0])
                    strands.append(aln_block[sp][3])

                if aln_block[sp][5][pos] == '-':
                    gap_count[sp] += 1
                    if indel is False:
                        positions.append('NA')
                else:
                    if indel is False:
                        # catch revcomp position here
                        spp_strand = aln_block[sp][3]

                        current_maf_pos = int(aln_block[sp][1]) + pos - gap_count[sp]
                        if spp_strand == '+':
                            positions.append(str(current_maf_pos))
                            if sp == ref:
                                ref_start = current_maf_pos

                        # adjust current relative block pos to genomic pos if on negative strand
                        else:
                            seq_len = int(aln_block[sp][4])
                            rev_coord_fix = seq_len - current_maf_pos - 1
                            positions.append(str(rev_coord_fix))

            else:
                sites[sp] += '?'
                if indel is False:
                    chroms.append('?')
                    strands.append('?')
                    positions.append('?')

    # print final record in block
    if start_print is True:  # deals with instances where whole block is an indel
        ref_end = ref_start + (len(sites[ref]) - sites[ref].count('-'))

        bed_line_str = [str(s) for s in [bed_line[0], ref_start, ref_end, bed_line[3]]]

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


def main():

    parser = argparse.ArgumentParser(description="Convert a whole genome alignment in MAF format to BED format "
                                                 "following the coordinates of one of the species in the alignment ")
    parser.add_argument('-i', '--infile',
                        dest='infile',
                        required=True,
                        help="Whole genome alignment in MAF format (Compressed)")
    parser.add_argument('-r', '--ref_species',
                        dest='ref_species',
                        required=True,
                        help="Name of reference species (as it appears in the MAF file")
    parser.add_argument('-c', '--chromosome',
                        dest='ref_chrom',
                        required=True,
                        help="Specify which chromosome to extract. This script only extracts one to BED format one "
                             "chromosome in each run")

    args = parser.parse_args()

    ref_species = args.ref_species

    # generate list of species in maf file, reference first, then alphabetical
    spp_grep = "zgrep ^s " + args.infile + " | head -n 50000 | cut -d '.' -f 1 | cut -d ' ' -f 2 | less -S | sort -u"
    species_list = subprocess.Popen(spp_grep, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    species_list.remove(ref_species)
    species_list = [ref_species] + sorted(species_list)

    ref_chrom = args.ref_chrom

    with gzip.open(args.infile, 'r') as infile:
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
                align_block[species[1].split('.')[0]] = ['.'.join(species[1].split('.')[1:])] + species[2:]

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
