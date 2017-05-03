#!/usr/bin/env python

from __future__ import print_function
import sys


def main():
    for line in sys.stdin:
        seqs = line.split('\t')[7].split(',')
        is_masked = set([x.islower() for x in seqs])

        if len(is_masked) == 1 and True in is_masked:
            print(line)

if __name__ == '__main__':
    main()
