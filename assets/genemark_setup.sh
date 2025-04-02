#!/bin/bash
tool_link=$1
key_link=$2

# download interproscan and genemark
wget $tool_link

# untar the directory
tar -xvzf gmes_linux_64_4.tar.gz

# rename it
mv gmes_linux_64_4 gmes_linux_64
wget $key_link
gunzip gm_key_64.gz

# copy key
cp gm_key_64 gmes_linux_64/.gm_key

# in genemark, run script to fix perl paths inside container
cd gmes_linux_64
perl change_path_in_perl_scripts.pl "/opt/conda/envs/Annotate_NLR/bin/perl"
