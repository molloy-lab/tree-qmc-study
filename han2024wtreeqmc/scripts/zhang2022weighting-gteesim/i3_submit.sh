#!/bin/bash

exit

# STUDY PARAMETERS
REPLS=( $(seq -f "%02g" 1 50) )   # Replicates

for REPL in ${REPLS[@]}; do
MODL="model.$REPL"
echo "Submitting $MODL ..."

NAME="i3.$MODL"

sbatch \
    --job-name="$NAME" \
    --output="$NAME.%j.out" \
    --error="$NAME.%j.err" \
    --export=REPL="$REPL" \
i3_drive.sbatch
done

