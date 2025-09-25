#!/usr/bin/env python3

"""
Given an input tsv of input gff files or anticipated gff files, parse
the input:
    - [x] confirm helixer, stringtie, transdecoder, miniprot inputs are
              labelled correctly with hx, st, tr, and mp
    - [x] if input gff provided for the above expected fields, check
              that the path(s) exists locally
    - [x] if no input gffs in hx, str, tr, and mp, fill these fields
              with the generic names the pipeline will produce
    - [x] additional gff files exist on the system
    - [x] names are all unique
    - [x] if st or hx is skipped on run, make sure the gffs are present.
    - [x] modifiers are empty or numeric

Usage:
python3 parse_input.py <mikado design> <skip_st> <skip_hx> <nlrs> <lo_genome> <lo_gff>

Mikado configuration file
https://mikado.readthedocs.io/en/stable/Tutorial/

1. filepath
2. id (name used in gff output)
3. strand-specific?
4. modifier (can be positive or negative, missing means no modification)
5. reference?
6. exclude redundant models?
7. strip CDS of faulty models? (False will discard them)
8. skip chimera split routine? (False will proceed normally)

Example for ragnarok base pipeline (using generic gff file names):
helixer.gff3    hx    True        False    False
stringtie.gtf    st    True    1    False    True
transcripts.fasta.transdecoder.genome.gff3    tr    False    -0.5    False    False
aa_miniprot.gff    mp    True    1    False    False
"""

import os
import numpy as np
import pandas as pd
import sys


def get_input():
    """
    Load the tsv input and apply column labels based on input length.
    """
    colnames = ["file",
                "id",
                "strand",
                "modifier",
                "reference",
                "exclude_redundant",
                "strip_cds",
                "skip_chimera"]

    df = pd.read_csv(sys.argv[1], sep="\t", header=None)
    colnames = colnames[:df.shape[1]]
    df.columns = colnames

    return df, colnames


def expected_ids(df):
    """
    Assert that the id field contains the expected input ids.
    The helixer and stringtie, transdecoder, miniprot steps are by
    design included in the annotation pipeline.
    """
    exp_ls = ["hx", "st", "tr", "mp"]
    for i in exp_ls:
        assert df["id"].str.contains(i).any(), f"{i} not found"


def unique_ids(df):
    """
    Assert that the id field contains only unique ids.
    """
    id_ls = df["id"].to_list()
    uniq_id_ls = df["id"].unique().tolist()
    assert len(id_ls) == len(uniq_id_ls), f"{id_ls} duplicates found"


def fields_okay(df, colnames):
    """
    Confirm that input fields are bool or numeric as needed.
    If NaN present, replace with empty 0 for assert, then empty string.
    """
    # df["modifier"] = df["modifier"].replace(np.NaN, "0")
    df["modifier"] = df["modifier"].replace(np.nan, "0")

    for i in colnames[2:]:
        if i == "modifier":
            assert pd.to_numeric(df["modifier"], errors='coerce').notna().all(), f"modifier must be numeric"
        else:
            assert df[i].isin([True, False]).all(), f"strand must be bool for {i}"
    df["modifier"] = df["modifier"].replace("0", "")

    return df


def expected_gffs(df):
    """
    Assert that the ragnarok input gffs contain existing files, if user
    has input them.
    If hx is present, check for filepath. If missing, store in hx_ls.
    If st, tr, or mp present, assert all three have filepaths or fail as
    these three files are progressively built upon star output (store
    missing files in st_ls).
    Returns if stringtie is all present, helixer is present.
    """
    exp_ls = ["hx", "st", "tr", "mp"]
    st_ls = []
    hx_ls = []
    for i in exp_ls:
        gff = df.loc[df["id"] == i, "file"].squeeze()
        if pd.isnull(gff):
            if i in ["st", "tr", "mp"]:
                st_ls.append(i)
            else:
                hx_ls.append(i)
        else:
            assert os.path.exists(gff), f"{gff} path not found"

    assert len(st_ls) in {0, 3}, f"st, tr, mp must be present"

    return len(st_ls) == 0, len(hx_ls) == 0


def confirm_st_hx_status(st_present, hx_present):
    """
    If stringtie gffs detected, confirm this matches user input for the
    skip_st param. If helixer gff detected, confirm this matches the
    user input for skip_hx.
    """
    skip_st = True if sys.argv[2] == "true" else False
    skip_hx = True if sys.argv[3] == "true" else False
    assert st_present == skip_st, f"skip_st is {skip_st} but gff presence is {st_present}"
    assert hx_present == skip_hx, f"skip_hx is {skip_hx} but gff presence is {hx_present}"


def confirm_lo_status(df):
    """
    If user has lo_genome defined, assert that no input gff for lo id is
    present.
    """
    lo_genome = False if sys.argv[5] == "false" else sys.argv[5]
    lo_gff = False if sys.argv[6] == "false" else sys.argv[6]
    print(lo_genome)
    print(lo_gff)
    lo_id = (df["id"] == "lo").any()
    lo_present = df.loc[df["id"] == "lo", "file"].notna().any()

    if lo_id:
        assert lo_present != True, f"lo id found but reserved for use with params.lo_genome params.lo_gff, please rename id"
        assert os.path.exists(lo_genome), f"params.genome input but path to fasta not found"
        assert os.path.exists(lo_gff), f"params.gff input but path to gff not found"
    else:
        assert lo_genome is False
        assert lo_gff is False

    return lo_present


def confirm_nlr_status(df):
    """
    If user has nlrs == true, assert that no input gff for nlr id is
    present.
    """
    skip_nlr = False if sys.argv[4] == "true" else True
    nlr_id = (df["id"] == "nlr").any()
    nlr_present = df.loc[df["id"] == "nlr", "file"].notna().any()

    if nlr_id:
        assert skip_nlr != True, f"nlr id found and should be empty when params.nlr == true"
        assert nlr_present != True, f"nlr id found but reserved for use with params.nlr, please rename id"

    if not skip_nlr:
        assert nlr_id == True, f"nlr id needs to be defined in params.design"

    return nlr_present


def confirm_hx_um_status(df):
    """
    If user has mask_helixer == consensus, assert that no input gff for
    hx_um id is present.
    """
    skip_hx_um = False if sys.argv[7] == "consensus" else True
    hx_um_id = (df["id"] == "hx_um").any()
    hx_um_present = df.loc[df["id"] == "hx_um", "file"].notna().any()

    if hx_um_id:
        assert skip_hx_um != True, f"hx_um id found and should be empty when params.mask_helixer == consensus"
        assert hx_um_present != True, f"hx_um id found but reserved for use with params.mask_helixer, please rename id"

    if not skip_hx_um:
        assert hx_um_id == True, f"hx_um id needs to be defined in params.design"

    return hx_um_present


def write_existing(df):
    existing_files = df["file"].to_list()

    with open("gff_paths.csv", "w") as o:
        for i in existing_files:
            if pd.isnull(i):
                pass
            else:
                o.write(f"{i}\n")


def adjust_names(df, st_present, hx_present, nlr_present, lo_present, hx_um_present):
    """
    Once filepaths have been confirmed to exist or not, a table of all
    files as they will appear once symlinked to mikado2 process is
    created as the configuration file.
    """
    if not st_present:
        df.loc[df["id"] == "st", "file"] = "stringtie.gff"
        df.loc[df["id"] == "tr", "file"] = "transcripts.fasta.transdecoder.genome.gff3"
        df.loc[df["id"] == "mp", "file"] = "aa_miniprot.gff"
    if not hx_present:
        df.loc[df["id"] == "hx", "file"] = "helixer.gff3"
    if not hx_um_present:
        df.loc[df["id"] == "hx_um", "file"] = "helixer_um.gff3"
    if not lo_present:
        df.loc[df["id"] == "lo", "file"] = "liftoff.gff3"
    if not nlr_present:
        df.loc[df["id"] == "nlr", "file"] = "fpnlr_braker_aa_NLR.gff3"

    df["file"] = df["file"].apply(os.path.basename)

    return df


def main():
    df, colnames = get_input()
    expected_ids(df)
    unique_ids(df)
    df = fields_okay(df, colnames)
    st_present, hx_present = expected_gffs(df)
    confirm_st_hx_status(st_present, hx_present)
    lo_present = confirm_lo_status(df)
    nlr_present = confirm_nlr_status(df)
    hx_um_present = confirm_hx_um_status(df)
    write_existing(df)
    df = adjust_names(df, st_present, hx_present, nlr_present, lo_present, hx_um_present)
    df.to_csv("mikado.tsv", sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
