#!/usr/bin/env python3

"""
Given an input tsv of input gff files or anticipated gff files, parse
the input:
	- [x] confirm helixer, stringtie, transdecoder, miniprot inputs are
	  labelled correctly with hx, st, tr, and mp
	- [ ] if input gff provided for the above expected fields, check
	  that the path(s) exists locally
	- [ ] if no input gffs in hx, str, tr, and mp, fill these fields
		  with the generic names the pipeline will produce
	- [ ] additional gff files exist on the system
	- [ ] names are all unique
	- [ ] modifiers are empty or numeric


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
helixer.gff3	hx	True		False	False
10kIntron_stringtie.gtf	st	True	1	False	True
transcripts.fasta.transdecoder.genome.gff3	tr	False	-0.5	False	False
aa_miniprot.gff	mp	True	1	False	False
"""

import os
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
	df.columns = colnames[:df.shape[1]]

	return df


def expected_ids(df):
	"""
	Assert that the id field contains the expected input ids.
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


def adjust_names(df, st_present, hx_present):
	gff_ls = df["id"].unique().tolist()
	if not st_present:
		df.loc[df["id"] == "st", "file"] = "10kIntron_stringtie.gtf"
		df.loc[df["id"] == "tr", "file"] = "transcripts.fasta.transdecoder.genome.gff3"
		df.loc[df["id"] == "mp", "file"] = "aa_miniprot.gff"
	if not hx_present:
		df.loc[df["id"] == "hx", "file"] = "helixer.gff3"
	df = df["file"].apply(os.path.basename)

	return df


def main():
	df = get_input()
	expected_ids(df)
	unique_ids(df)
	st_present, hx_present = expected_gffs(df)
	print(st_present)
	print(hx_present)
	df = adjust_names(df, st_present, hx_present)
	print(df)


if __name__ == "__main__":
	main()
