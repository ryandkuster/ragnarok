#!/usr/bin/env python3

"""
Open the input protein fasta file (stdin1). If encountering a non-IUPAC
residue, replace with an X (as these are not alignments).
Write out to stdin2 file.
"""

import math
import sys


braker_prots = "AaRrNnDdCcEeQqGgHhIiLlKkMmFfPpSsTtWwYyVvBbZzJjOoUuXx"


def main():
    """
    Open the fasta protein file, write lines if they are headers, else
    send to check_prots functions.
    """
    with open(sys.argv[1]) as f, open(sys.argv[2], "w") as o:
        for line in f:
            if line.startswith('>'):
                o.write(line)
            else :
                line = line.rstrip()
                line = check_prot(line)
                o.write(line)


def check_prot(line):
    """
    Check the amino acids and replace non "braker_prots" residues with X.
    """
    new_line = ''
    for aa in line:
        if aa in braker_prots:
            new_line += aa
        else:
            new_line += "X"
    return f"{new_line}\n"


if __name__ == "__main__":
    main()
