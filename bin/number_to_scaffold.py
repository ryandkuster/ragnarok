#!/usr/bin/env python3

"""
Corresponding script to run using output from scaffold_to_number.py.

Rename fasta using original key:value relationship of headers.
Open the input reference fasta file (stdin1).
Write out renamed fasta to stdin2.
Open the input dictionary to rename with (stdin3).
"""

import json
import sys


def main():
    """
    Process fasta file and rename scaffolds (headers) to the
    corresponding value in the key_dt.
    """
    in_fasta  = sys.argv[1]
    out_fasta = sys.argv[2]
    key_file  = sys.argv[3]

    with open(key_file, "r") as f:
        key_dt = json.load(f)

    with open(in_fasta) as f, open(out_fasta, "w") as o:
        for line in f:
            if line.startswith(">"):
                # number = int(line.rstrip()[1:])
                number = line.rstrip()[1:]
                scaffold = key_dt[number]
                o.write(f">{scaffold}\n")
            else:
                o.write(line)


if __name__ == "__main__":
    main()
