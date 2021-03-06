#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


def species_in_block(positon_column):

    """
    takes the coordinates column of wga.bed and returns number of spp covered
    :param positon_column: list
    :return: int
    """

    no_spp_align = len(positon_column)
    spp_missing = positon_column.count('?')
    spp_in_block = no_spp_align - spp_missing
    return spp_in_block


def unique_to_ref(sequences):

    """
    takes the sequence column of a wga.bed and returns true if INDEL is reference specific
    :param sequences: list
    :return: bool
    """

    ref = sequences[0]
    outgroups = sequences[1:]

    # check for seq conservation in outgroup
    for spp in outgroups:
        if spp.count('-') != outgroups[0].count('-'):
            return False

    # check ref not same as outgroups
    if ref.count('-') == outgroups[0].count('-'):
        return False
    else:
        return True


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-max_length', help='Maximum INDEL length to extract', default=None, type=int)
    parser.add_argument('-min_coverage', help='Minimum species coverage', default=1, type=int)
    parser.add_argument('-ref_specific', help='Restricts output to INDELs only found in reference species, and '
                                              'conserved in non reference species', default=False, action='store_true')
    parser.add_argument('-lengths', help='Comma separated list of lengths to extract', default=None)
    args = parser.parse_args()

    # variables
    max_len = args.max_length
    min_spp = args.min_coverage
    ref_specif = args.ref_specific

    # extraction
    for orig_line in sys.stdin:
        line = orig_line.rstrip('\n').split('\t')
        chromo, start, end, strand, spp, all_chromo, all_pos, sites, score = line[0], line[1], line[2], line[3], \
            line[4].split(','), line[5].split(','), line[6].split(','), line[7].split(','), line[8]

        # skips non INDELs
        if '-' not in ''.join(sites):
            continue

        else:
            # length checks
            if len(sites[0])-1 > max_len and max_len is not None:
                continue
            if args.lengths is not None:
                lengths = [int(x) for x in args.lengths.split(',')]
                if len(sites[0])-1 not in lengths:
                    continue

            # coverage check
            if species_in_block(all_pos) < min_spp:
                continue

            # ref specific?
            if ref_specif is True and unique_to_ref(sites) is False:
                continue

            print(orig_line.rstrip('\n'))

if __name__ == '__main__':
    main()
