 #!/usr/bin/env bash

set -x
set -e

script_dir="$(dirname $0)"

nextflow run "$script_dir"/.. \
    -stub -resume -profile gh \
    --sample_sheet "$script_dir"/sample-sheet.csv \
    --inputs "$script_dir"/inputs \
    --outputs "$script_dir"/outputs
