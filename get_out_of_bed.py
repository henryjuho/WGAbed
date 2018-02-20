#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam


def get_start_index(start, seq_list, pos_list):

    """
    finds where desired sequence starts in interval
    :param start: int
    :param seq_list: list
    :param pos_list: list
    :return: int
    """

    # if start or region is the start of the bed line
    if start == pos_list[0]:
        return 0

    # if this is a continuation of the region
    elif start < pos_list[0]:
        return 0

    # if the region starts within the bed line
    else:
        extra_len = start - pos_list[0]

        # find pos in ref seq
        base_count = 0
        char_count = 0
        for b in seq_list[0]:
            char_count += 1
            if b != '-':
                base_count += 1
                if base_count == extra_len:
                    break

        return char_count


def get_end_index(stop, seq_list, pos_list):

    """
    finds where desired sequence ends in interval
    :param stop: int
    :param seq_list: list
    :param pos_list: list
    :return: int
    """

    # if end of regions is end of bed interval
    if stop == pos_list[1]:
        return len(seq_list[0])

    # if end of region is after bed interva;
    elif stop > pos_list[1]:
        return len(seq_list[0])

    # if end of region is within feature
    else:
        extra_len = pos_list[1] - stop

        # find pos in ref seq
        base_count = 0
        char_count = 0
        for b in seq_list[0][::-1]:
            char_count += 1
            if b != '-':
                base_count += 1
                if base_count == extra_len:
                    break

        return -char_count


def merge_sequences(seq_list):

    """
    merges a list of sequence lists
    :param seq_list: list
    :return: list
    """
    merged_align = [''.join(y) for y in zip(*seq_list)]
    return merged_align


def rm_ins_rel_ref(spp, seq):

    """
    strips out insertions relative to the ref
    :param spp: tuple
    :param seq: list
    :return: list
    """

    n_spp = len(spp)

    if n_spp == 0:
        return []

    trimmed_seqs = ['' for x in spp]

    for pos in range(0, len(seq[0])):

        ref_base = seq[0][pos]

        if ref_base == '-':
            continue

        spp_seqs = [y[pos] for y in seq]

        trimmed_seqs = [trimmed_seqs[i] + spp_seqs[i] for i in range(0, n_spp)]

    return trimmed_seqs


def intersect2align(chromo, start, end, wga_bed, ins_rel_ref=True):

    """
    get subregions from wga bed
    :param chromo: str
    :param start: int
    :param end: int
    :param wga_bed: str
    :param ins_rel_ref: bool
    :return: tuple, list
    """
    var_align = [x for x in wga_bed.fetch(chromo, start, end, parser=pysam.asTuple())]

    counter = 0
    concat_seqs = []
    previous_end = 0
    seq_ids = ()
    for bed_row in var_align:
        counter += 1
        positions = [int(x) for x in bed_row[1:3]]
        seq_ids = tuple(bed_row[4].split(','))
        sequences = bed_row[7].split(',')
        no_species = len(sequences)

        # work out gap adjustment
        if counter != 1 and previous_end != positions[0]:
            gap_fill = [''.join(['N' for i in range(previous_end, positions[0])]) for j in range(0, no_species)]
        else:
            gap_fill = ['' for j in range(0, no_species)]

        # if first bed interval returned
        if counter == 1:
            concat_seqs = ['' for j in range(0, no_species)]

            # if beginning of region not in alignment
            if positions[0] > start:
                missing_len = positions[0] - start
                seq_correction = [''.join(['N' for i in range(0, missing_len)]) for j in range(0, no_species)]
                concat_seqs = seq_correction

        # get section of bed interval needed
        start_index = get_start_index(start, sequences, positions)
        stop_index = get_end_index(end, sequences, positions)
        interval_seqs = [x[start_index:stop_index] for x in sequences]
        concat_seqs = merge_sequences([concat_seqs, gap_fill, interval_seqs])

        # if last bed interval returned
        if counter == len(var_align):
            # if end of region not in alignment
            if positions[1] < end:
                missing_len = end - positions[1]
                seq_correction = [''.join(['N' for i in range(0, missing_len)]) for j in range(0, no_species)]
                concat_seqs = merge_sequences([concat_seqs, gap_fill, sequences, seq_correction])

        previous_end = positions[1]

    if not ins_rel_ref:
        concat_seqs = rm_ins_rel_ref(seq_ids, concat_seqs)

    return seq_ids, concat_seqs


def main():
    # arguments
    parser = argparse.ArgumentParser(description='Utility to convert subregions of an whole genome alignment '
                                                 'bed file into different formats, eg) fasta')
    parser.add_argument('-wb', '--wga_bed', help='Whole genome alignment bedfile', required=True)
    parser.add_argument('-q', '--query', help='Extraction coordinates in form chr:start-end', required=True)
    parser.add_argument('-f', '--format', help='Format to output', required=True, choices=['fasta', 'phylip'])
    args = parser.parse_args()

    # variables
    wb = pysam.TabixFile(args.wga_bed)
    q = (args.query.split(':')[0], int(args.query.split('-')[0].split(':')[1]), int(args.query.split('-')[1]))
    out_format = args.format

    # get data out
    extracted_data = intersect2align(q[0], q[1], q[2], wb)

    if out_format == 'fasta':
        for i in range(0, len(extracted_data[0])):
            print('>' + extracted_data[0][i])
            for q in range(0, len(extracted_data[1][i]), 60):
                print(extracted_data[1][i][q: q+60])
    if out_format == 'phylip':
        print('\t{}\t{}'.format(len(extracted_data[0]), len(extracted_data[1][0])))
        for i in range(0, len(extracted_data[0])):
            print(extracted_data[0][i])
            for q in range(0, len(extracted_data[1][i]), 60):
                print(extracted_data[1][i][q: q+60])

if __name__ == '__main__':
    main()
