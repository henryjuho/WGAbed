#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip
import sys


def spp_position(spp_col, target_spp):

    """
    takes the spp column string from wga.bed and returns the index of the query species
    :param spp_col: str
    :param target_spp: str
    :return: int
    """
    spp_order = spp_col.split(',')
    try:
        return spp_order.index(target_spp)
    except ValueError:
        sys.exit('Query species is not in wga.bed, check your spelling!')


def wgasite_in_query(wga_start, wga_end, query_set):

    """
    takes the start and end coordinates of a non ref site and identifies if it is in the query bed file
    :param wga_start: int
    :param wga_end: int
    :param query_set: set
    :return: bool
    """
    wga_positions = range(wga_start, wga_end)
    for site in wga_positions:
        if site not in query_set:
            return False
    return True


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed_non_ref',
                        help='.bed file with non refernce spp coordinates to intersect with wga.bed file',
                        required=True)
    parser.add_argument('-q', '--query_species',
                        help='species that the query bed file has coordinates for',
                        required=True)
    parser.add_argument('-c', '--chromosome',
                        help='target chromosome',
                        required=True)
    args = parser.parse_args()

    # variables
    query_spp = args.query_species
    query_chromo = args.chromosome
    non_ref_coords = set(sum([range(int(x.split()[1]), int(x.split()[2])) for x in gzip.open(args.bed_non_ref) if
                              x.split()[0] == query_chromo], []))

    # loop through wga.bed lines from stdin
    for orig_line in sys.stdin:
        line = orig_line.rstrip('\n').split('\t')
        spp_index = spp_position(line[4], query_spp)
        chromo, pos, sites = line[5].split(',')[spp_index], line[6].split(',')[spp_index], line[7].split(',')[spp_index]

        # skip sites where query spp is not covered
        if chromo == '?' or pos == '?' or sites == '?':
            continue

        start = int(pos)
        end = start + (len(sites) - sites.count('-'))

        # skip non target chromos
        if chromo != query_chromo:
            continue

        # see if site is in query skip if not
        if wgasite_in_query(start, end, non_ref_coords) is False:
            continue

        # print all passing lines
        print(orig_line.rstrip('\n'))

if __name__ == '__main__':
    main()
