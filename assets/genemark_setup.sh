#!/bin/bash

# download interproscan and genemark
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_aqg0O/gmes_linux_64_4.tar.gz

# untar the directory
tar -xvzf gmes_linux_64_4.tar.gz

# rename it
mv gmes_linux_64_4 gmes_linux_64
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_aqg0O/gm_key_64.gz
gunzip gm_key_64.gz

# copy key
cp gm_key_64 gmes_linux_64/.gm_key

# in genemark, run script to fix perl paths inside container
cd gmes_linux_64
perl change_path_in_perl_scripts.pl "/opt/conda/envs/Annotate_NLR/bin/perl"
