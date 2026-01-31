#!/bin/bash

MAIN_DIR="$1"

echo -e "File\tRead_Count"

find "$MAIN_DIR" -type f -name "*.bam" | while read -r BAM; do
    SAMPLE=$(basename "$BAM")
    READ_COUNT=$(samtools view -c -F 260 "$BAM")
    echo -e "${SAMPLE}\t${READ_COUNT}"
done