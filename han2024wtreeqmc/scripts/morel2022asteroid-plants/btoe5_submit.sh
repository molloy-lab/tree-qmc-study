#!/bin/bash

exit

INFO="plants"

sbatch \
    --job-name="btoe5.$INFO" \
    --output="btoe5.$INFO.%j.out" \
    --error="btoe5.$INFO.%j.err" \
btoe5_drive.sbatch

