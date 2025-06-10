#!/usr/bin/env bash

# If for some reason your server can't automatically pull apptainer
# images, provide the same scratch directory you plan to use for all
# ragnarok runs and this script should download everything.

SCRATCHDIR=$1
sing_dir=$SCRATCHDIR/singularity
mkdir -p $sing_dir

script_dir=$( realpath $( dirname $0 ) )
rag_dir="${script_dir%/*}"

#grep "container" $rag_dir/conf/containers.config | awk '{print $3}'
for image in $( grep "container" $rag_dir/conf/containers.config | awk '{print $3}' ) ; do 
    image=${image#\"}
    image=${image%\"}
    if [[ $image == docker* ]]; then
        echo $image
    else
        image="docker://${image}"
        echo $image
    fi

    sif_name=$( basename $image )
    sif_name=${sif_name}.sif

    apptainer pull --dir $sing_dir $sif_name $image
done

echo -e "Downloaded images to ${SCRATCHDIR}/singularity.\n"
echo -e "Use SCRATCHDIR=$SCRATCHDIR in nextflow runs.\n"
