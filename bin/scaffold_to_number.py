#!/usr/bin/env python3

"""
Open the input reference fasta file (stdin1).
Check if input is gzipped and handle.
For each header, rename in order of appearance to integer.
Write out renamed fasta to stdin2.
Save original key:value relationship of old:new names to file (stdin3)
"""

import gzip
import json
import sys


def gzip_test(fasta):
    """
    Check if the input fasta is gzipped and return True if it is.
    """
    try:
        with open(fasta) as f:
            f.readline()
        compressed = False
    except UnicodeDecodeError:
        compressed = True
    return compressed


def rename_scaffolds(f, out_fasta):
    """
    Process fasta file and rename scaffolds (headers) to current value
    of counter. Each time, store counter and original header to key_dt.
    Return key_dt.
    """
    counter = 0
    key_dt = {}

    with open(out_fasta, "w") as o:
        for line in f:
            if line.startswith(">"):
                counter += 1
                o.write(f">{counter}\n")
                key_dt[counter] = line.rstrip()[1:]
            else:
                o.write(line)
    return key_dt



def main():
    in_fasta  = sys.argv[1]
    out_fasta = sys.argv[2]
    key_file  = sys.argv[3]

    if gzip_test(in_fasta):
        with gzip.open(in_fasta, "rt") as f:
            key_dt = rename_scaffolds(f, out_fasta)
    else:
        with open(in_fasta) as f:
            key_dt = rename_scaffolds(f, out_fasta)

    with open(key_file, "w") as f:
        json.dump(key_dt, f)


if __name__ == "__main__":
    main()
