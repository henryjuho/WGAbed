#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--species', help='Species coordinates to use, currently reference only', required=True)
    parser.add_argument('-r', '--region', help='Region to extract in form chr:start-stop', required=True)
    args = parser.parse_args()

    # variables
    spp = args.species
    region = args.region.replace('-', ':').split(':')
    region = [region[0], int(region[1]), int(region[2])]
    # process piped maf
    align_block = []
    printable = False

    for line in sys.stdin:
        # print header
        if line.startswith('#'):
            print(line.rstrip())

        # write old block and store new
        elif line.startswith('a'):
            if printable:
                print(''.join(align_block))
            align_block = [line]
            printable = False
        elif line.startswith('s'):
            species = line.split()[1].split('.')[0]
            chromo = line.split()[1].split('.')[1]

            # if reference spp and chromo
            if spp == species and chromo == region[0]:
                pos = int(line.split()[2])
                if region[1] <= pos < region[2]:
                    printable = True
            align_block.append(line)

    # write last block if in region
    if printable:
        print(''.join(align_block))

if __name__ == '__main__':
    main()
