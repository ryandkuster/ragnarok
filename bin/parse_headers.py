#!/usr/bin/env python3

"""
Given a masked and unmasked genome file, confirm that the names of all headers
are identical between the two.

Usage:
python3 parse_headers.py <genome> <genome_masked>

1. filepath unmasked genome
2. filepath masked genome
"""

import gzip
import os
import sys


def get_headers(fna):
    headers = []
    try:
        with open(fna) as f:
            for line in f:
                if line.startswith(">"):
                    headers.append(line.rstrip()[1:])
    except UnicodeDecodeError:
        with gzip.open(fna, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    headers.append(line.rstrip()[1:])
    headers.sort()
    return headers


def compare_headers(headers_1, headers_2):
    for i, j in zip(headers_1, headers_2):
        assert i == j, "sorted headers do not match between genomes"
    print("all headers passed")


def main():
    fna_1 = sys.argv[1]
    fna_2 = sys.argv[2]
    headers_1 = get_headers(fna_1)
    headers_2 = get_headers(fna_2)
    compare_headers(headers_1, headers_2)


if __name__ == "__main__":
    main()

