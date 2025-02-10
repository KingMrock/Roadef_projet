#!/bin/bash
# run_all_instances.sh
# Time limit for each instance (30 minutes = 1800 seconds)
TIME_LIMIT=1800

# Base directory containing the instance files.
BASE_DIR="/data/storage/PRJ"

# Path to the solver (adjust if needed).
SOLVER_SCRIPT="solver.py"

# Find all JSON files recursively under BASE_DIR.
find "$BASE_DIR" -type f -name "*.json" | while read INSTANCE_FILE; do
    # Extract the base name (without extension) from the instance file.
    INSTANCE_BASENAME=$(basename "$INSTANCE_FILE" .json)
    # Define the solution file name (same folder as the instance).
    SOLUTION_FILE="$(dirname "$INSTANCE_FILE")/${INSTANCE_BASENAME}.txt"
    
    echo "-----------------------------------------"
    echo "Running solver for instance: $INSTANCE_FILE"
    echo "Solution file: $SOLUTION_FILE"
    
    # Launch the solver with the specified time limit.
    python3 "$SOLVER_SCRIPT" "$INSTANCE_FILE" -s "$SOLUTION_FILE" --time "$TIME_LIMIT"
    
    echo "Finished instance: $INSTANCE_FILE"
    echo "-----------------------------------------"
done
