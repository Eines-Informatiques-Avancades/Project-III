#!/bin/bash

for file in thermodynamics_*.dat; do
    echo "Processing file: $file"
    export FILE="$file"

    # Your processing logic goes here
    # You can use $FILENAME in your processing steps
    python3 compute_average.py >> results.dat
done

