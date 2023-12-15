#!/bin/bash

set -euxo pipefail

snakemake --printshellcmds --use-conda --conda-prefix ../conda-env --cores 40 --rerun-triggers mtime --default-resources mem_mb=172032 --set-resources panvc3_generate_founder_sequences:mem_mb=81920 --dry-run -- all
