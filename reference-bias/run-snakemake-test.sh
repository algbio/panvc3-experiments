#!/bin/bash

set -euxo pipefail

snakemake --printshellcmds --use-conda --conda-prefix ../conda-env --cores 40 --rerun-triggers mtime --dry-run -- all
