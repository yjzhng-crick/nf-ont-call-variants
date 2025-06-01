 #!/usr/bin/env bash

set -x
set -e

script_dir="$(dirname $0)"
DOCKER=${1:-no}

if [ "$DOCKER" == "docker" ]
then
    export NXF_CONTAINER_ENGINE=docker
    docker_flag='-stub -with-docker -profile gh'
else
    export SINGULARITY_FAKEROOT=1
    export SINGULARITY_TMPDIR="$script_dir/singularity-tmp"
    docker_flag=''
fi

mkdir -p "$SINGULARITY_TMPDIR"
nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --sample_sheet "$script_dir"/sample-sheet.csv \
    --inputs "$script_dir"/inputs \
    --outputs "$script_dir"/outputs
