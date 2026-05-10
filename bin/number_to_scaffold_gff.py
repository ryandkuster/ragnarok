#!/usr/bin/env python3

"""
Corresponding script to run using output from scaffold_to_number.py.

Rename GFF3 seqnames using the original key:value relationship of headers.
Takes an input GFF3 file (--gff), an output path (--out_gff), and a JSON
key file (--json) mapping numeric seqnames to scaffold identifiers.
"""

import argparse
import json
import sys


def main():
    """
    Process gff file and rename scaffolds (headers) to the
    corresponding value in the key_dt.
    """
    parser = argparse.ArgumentParser(
        description="Rename GFF3 seqnames using a number-to-scaffold JSON key file."
    )
    parser.add_argument("--gff", required=True, help="Input GFF3 file.")
    parser.add_argument("--out_gff", required=True, help="Output GFF3 file.")
    parser.add_argument("--json", required=True, help="JSON key file mapping numbers to scaffold names.")
    args = parser.parse_args()

    in_gff   = args.gff
    out_gff  = args.out_gff
    key_file = args.json

    with open(key_file, "r") as f:
        key_dt = {k: v.split()[0] for k, v in json.load(f).items()}

    with open(in_gff) as f, open(out_gff, "w") as o:
        for line in f:
            if not line.strip():
                continue
            if line.startswith("#"):
                o.write(line)
            else:
                fields = line.rstrip('\n').split('\t')
                number = fields[0]
                if number not in key_dt:
                    sys.exit(f"Error: seqname '{number}' not found in key file. Aborting.")
                fields[0] = key_dt[number]
                out_line = "\t".join(fields)
                o.write(f"{out_line}\n")


if __name__ == "__main__":
    main()
