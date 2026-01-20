#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    exit 1
fi

input_file="$1"
output_file="$2"

awk '
{
    n = NF
    # Separate all fields by a single space
    for (i=1; i<=n; i++) {
        printf "%s%s", $i, (i<n ? " " : "\n")
    }
}' "$input_file" > "$output_file"

echo "Processing complete. Cleaned tensor written to '$output_file'."
