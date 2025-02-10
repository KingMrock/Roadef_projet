#!/usr/bin/env python3
import os
import subprocess

# Set the time limit (in seconds) and base directory containing the JSON instances.
TIME_LIMIT = 1800  # 30 minutes
BASE_DIR = "/data/storage/PRJ"  # Adjust if needed
SOLVER_SCRIPT = "solver.py"      # Assumes solver.py is in the same folder as this script

def find_json_files(base_dir):
    """
    Recursively find and return a list of all .json files under base_dir.
    """
    json_files = []
    for root, dirs, files in os.walk(base_dir):
        for file in files:
            if file.endswith(".json"):
                json_files.append(os.path.join(root, file))
    return json_files

def run_solver(instance_file):
    """
    Run solver.py on a given instance file.
    The solution file will be in the same directory as the instance file,
    with the same base name and a .txt extension.
    """
    instance_dir = os.path.dirname(instance_file)
    instance_basename = os.path.splitext(os.path.basename(instance_file))[0]
    solution_file = os.path.join(instance_dir, instance_basename + ".txt")
    
    cmd = [
        "python3",
        SOLVER_SCRIPT,
        instance_file,
        "-s", solution_file,
        "--time", str(TIME_LIMIT)
    ]
    print("=========================================")
    print(f"Running solver for instance: {instance_file}")
    print(f"Solution file: {solution_file}")
    print("Command:", " ".join(cmd))
    
    subprocess.run(cmd)
    
    print(f"Finished instance: {instance_file}")
    print("=========================================")

def main():
    json_files = find_json_files(BASE_DIR)
    if not json_files:
        print(f"No JSON instance files found in {BASE_DIR}")
        return
    # Run the solver on each JSON file found.
    for instance_file in json_files:
        run_solver(instance_file)

if __name__ == "__main__":
    main()
