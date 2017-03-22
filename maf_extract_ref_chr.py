#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys


def main():

    parser = argparse.ArgumentParser(description="Extract a reference chromosome from a WGA MAF file")
    parser.add_argument('-c', '--chromosome',
                        required=True,
                        help="Specify which chromosome to extract")
    parser.add_argument('-H', '--header',
                        help='If specified will print header in output',
                        default=False,
                        action='store_true')
    args = parser.parse_args()

    chromo = args.chromosome
    header = args.header

    # loop through maf in stdin
    block = []
    spp_counter = 0
    ref_chr = 'none'

    for line in sys.stdin:

        # process header
        if line.startswith('#'):
            if header is True:
                print(line)
            else:
                continue

        # process block
        elif line.startswith('a'):
            block.append(line)

        elif line.startswith('s'):
            if spp_counter == 0:
                ref_chr = line.split()[1].split('.')[1]
            spp_counter += 1
            block.append(line)

        # print block
        else:
            if ref_chr == chromo:
                print(''.join(block))
            block = []
            spp_counter = 0

if __name__ == '__main__':
    main()
