#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse
from wga_bed_indels import species_in_block
from collections import OrderedDict


def out_group_agreement(sites):

    """
    returns true if outgroups agree
    :param sites: list(str, ...)
    :return: bool
    """

    out_groups = sites[1:]

    if len(set(out_groups)) == 1:
        return True

    else:
        return False


def main():

    # args
    parser = argparse.ArgumentParser()
    parser.add_argument('-callable',
                        help='Outputs number of sites that have all species covered and outgroup agreement',
                        default=False, action='store_true')
    args = parser.parse_args()

    out_dict = OrderedDict()

    # process file
    for orig_line in sys.stdin:
        line = orig_line.rstrip('\n').split('\t')
        chromo, start, end, strand, spp, all_chromo, all_pos, sites, score = line[0], line[1], line[2], line[3], \
            line[4].split(','), line[5].split(','), line[6].split(','), line[7].split(','), line[8]

        if chromo not in out_dict.keys():
            out_dict[chromo] = 0

        if args.callable:

            # coverage check
            if species_in_block(all_pos) < len(all_chromo):
                continue

            # ref specific?
            if out_group_agreement(sites) is False:
                continue

        out_dict[chromo] += 1

    for x in out_dict.items():
        print(x[0], x[1], sep='\t')

if __name__ == '__main__':
    main()
