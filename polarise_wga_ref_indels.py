#!/usr/bin/env python

from __future__ import print_function
import sys
from wga_bed_indels import unique_to_ref
import argparse


def indel_type(wga_line):

    """
    takes a wga bed line and returns indel type
    :param wga_line: str
    :return: bool
    """

    # get sequence column
    seqs = wga_line.split()[7].split(',')

    # catch non indel lines
    if len(seqs[0]) == 1:
        return 'non_indel'

    # process indel lines
    elif not unique_to_ref(seqs):
        return 'ambig_indel'

    # process ref specific indels
    else:
        ref_seq = seqs[0]

        if '-' in ref_seq:
            return 'deletion'

        else:
            return 'insertion'


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-indel_type', help='type of indel to output',
                        choices=['non_indel', 'ambig_indel', 'deletion', 'insertion'],
                        required=True)
    args = parser.parse_args()

    for line in sys.stdin:
        if indel_type(line) == args.indel_type:
            print(line.rstrip())


if __name__ == '__main__':
    main()
