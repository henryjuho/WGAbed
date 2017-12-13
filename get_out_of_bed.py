#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam


def merge_sequences(seq_list):

    """
    merges a list of sequence lists
    :param seq_list: list
    :return: list
    """
    merged_align = [''.join(y) for y in zip(*seq_list)]
    return merged_align


def intersect2align(chromo, start, end, wga_bed):

    """
    get subregions from wga bed
    :param chromo: str
    :param start: int
    :param end: int
    :param wga_bed: str
    :return:
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

        if counter == 1:
            # if beginning of region not in alignment
            if positions[0] > start:
                missing_len = positions[0] - start
                seq_correction = [''.join(['N' for i in range(0, missing_len)]) for j in range(0, no_species)]
                concat_seqs = merge_sequences([seq_correction, sequences])

            # if beginning of region is within first feature
            elif positions[0] < start:
                extra_len = start - positions[0]

                # find pos in ref seq
                base_count = 0
                char_count = 0
                for b in sequences[0]:
                    char_count += 1
                    if b != '-':
                        base_count += 1
                        if base_count == extra_len:
                            break

                concat_seqs = [x[char_count:] for x in sequences]

            else:
                concat_seqs = sequences

        elif counter == len(var_align):
            # if end of region not in alignment
            if positions[1] < end:
                missing_len = end - positions[1]
                seq_correction = [''.join(['N' for i in range(0, missing_len)]) for j in range(0, no_species)]
                concat_seqs = merge_sequences([concat_seqs, gap_fill, sequences, seq_correction])

            # if end of region is within last feature
            elif positions[1] > end:
                extra_len = positions[1] - end

                # find pos in ref seq
                base_count = 0
                char_count = 0
                for b in sequences[0][::-1]:
                    char_count += 1
                    if b != '-':
                        base_count += 1
                        if base_count == extra_len:
                            break

                concat_seqs = merge_sequences([concat_seqs, gap_fill, [x[:len(x)-char_count] for x in sequences]])

            else:
                concat_seqs = merge_sequences([concat_seqs, gap_fill, sequences])

        else:
            concat_seqs = merge_sequences([concat_seqs, gap_fill, sequences])

        previous_end = positions[1]

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
