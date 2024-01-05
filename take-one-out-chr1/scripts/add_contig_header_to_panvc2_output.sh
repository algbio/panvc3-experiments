#!/bin/bash

set -euxo pipefail

head -n 1 "${1}" > "${2}"
echo '##contig=<ID=1,length=249250621>' >> "${2}"
tail -n +2 "${1}" >> "${2}"
