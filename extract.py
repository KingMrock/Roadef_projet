#!/usr/bin/env python3
import os
import csv
import sys
import re

# Dictionary mapping instance IDs to the best objective value for 15 min.
best_known = {
    # A instances
    "A01": "1767,815611",
    "A02": "4671,376611",
    "A03": "848,1786111",
    "A04": "2085,876054",
    "A05": "635,2217857",
    "A06": "590,6235989",
    "A07": "2272,782274",
    "A08": "744,2932352",
    "A09": "1507,284784",
    "A10": "2994,848735",
    "A11": "495,2557702",
    "A12": "789,6349276",
    "A13": "1998,662162",
    "A14": "2264,124321",
    "A15": "2268,569150",
    # B instances
    "B01": "3986,202830",
    "B02": "4302,774528",
    "B03": "35279,530189",
    "B04": "34827,869811",
    "B05": "2397,099057",
    "B06": "4287,894340",
    "B07": "7564,034906",
    "B08": "7435,719048",
    "B09": "7491,753571",
    "B10": "10637,620000",
    "B11": "3626,271698",
    "B12": "37602,966667",
    "B13": "5024,492647",
    "B14": "11905,109948",
    "B15": "22566,003400",
    # C instances
    "C01": "8515,903773",
    "C02": "3541,653773",
    "C03": "33511,7",
    "C04": "37585,731132",
    "C05": "3166,889622",
    "C06": "8396,000943",
    "C07": "6083,270238",
    "C08": "11162,836",
    "C09": "5586,979245",
    "C10": "43342,488725",
    "C11": "5749,957352",
    "C12": "12721,133246",
    "C13": "42487,992826",
    "C14": "26467,221136",
    "C15": "39758,0275",
}

# Folder where checker outputs are stored
CHECKER_FOLDER = "checker_folder"

def extract_allowed_time(file_path):
    """Extract the allowed time from the file name.
    Expected filename format: L_NN_timeXXXX_log.txt, where XXXX is the allowed time.
    """
    filename = os.path.basename(file_path)
    match = re.search(r'time(\d+)_log', filename)
    return match.group(1) if match else ""

def extract_instance_id(file_path):
    """
    Extract the instance identifier from the file name.
    For checker file lookup, we expect the file to be named like A_01.txt.
    For example, from "A_01_time900_log.txt" or "A_01.json", we extract "A_01".
    """
    filename = os.path.basename(file_path)
    name, _ = os.path.splitext(filename)
    # Use a regex to capture the first part of the filename in the form "L_XX"
    m = re.match(r'([ABC]_\d{1,2})', name)
    if m:
        return m.group(1)
    return ""


def extract_checker_value(file_path, checker_folder=CHECKER_FOLDER):
    """
    Look for the corresponding checker output file in checker_folder.
    The checker file is expected to be named as instance_id.txt (e.g. A_01.txt).
    Extract the line starting with "Total objective" and return the value after the colon.
    """
    instance_id = extract_instance_id(file_path)
    checker_file = os.path.join(checker_folder, instance_id + ".txt")
    if os.path.exists(checker_file):
        with open(checker_file, "r") as f:
            for line in f:
                if line.lstrip().startswith("Total objective"):
                    parts = line.split(":", 1)
                    if len(parts) == 2:
                        return parts[1].strip()
    return ""

def extract_values(file_path):
    """Extracts values from a log file and augments with allowed time, best known objective,
    and the checker computed total objective value.
    """
    result = {
        "log_instance": "",
        "allowed_time": "",
        "time_to_build_model": "",
        "num_variables": "",
        "num_constraints": "",
        "total_solve_time": "",
        "time_to_first_solution": "",
        "first_solution_value": "",
        "gap_to_optimum": "",
        "best_solution_value": "",
        "best_obj_15min": "",
        "checker_total_objective": ""
    }
    
    # Save the full path or identifier from the file.
    result["log_instance"] = file_path
    result["allowed_time"] = extract_allowed_time(file_path)
    
    # Extract from the log file content.
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("Log for instance:"):
                result["log_instance"] = line.split(":", 1)[1].strip()
            elif line.startswith("Time to build model:"):
                result["time_to_build_model"] = line.split(":", 1)[1].strip().split()[0]
            elif line.startswith("Number of variables:"):
                result["num_variables"] = line.split(":", 1)[1].strip().split()[0]
            elif line.startswith("Number of constraints:"):
                result["num_constraints"] = line.split(":", 1)[1].strip().split()[0]
            elif line.startswith("Total solve time:"):
                result["total_solve_time"] = line.split(":", 1)[1].strip().split()[0]
            elif line.startswith("Time to first solution:"):
                result["time_to_first_solution"] = line.split(":", 1)[1].strip().split()[0]
            elif line.startswith("First solution value:"):
                result["first_solution_value"] = line.split(":", 1)[1].strip().split()[0]
            elif line.startswith("Gap to optimum:"):
                result["gap_to_optimum"] = line.split(":", 1)[1].strip().split()[0]
            elif line.startswith("Best solution value:"):
                result["best_solution_value"] = line.split(":", 1)[1].strip().split()[0]

    # Extract best known objective value for 15 minutes based on instance ID.
    instance_id = extract_instance_id(file_path)

    
    instance_id_lookup = instance_id.replace("_", "")
    if instance_id_lookup in best_known:
        result["best_obj_15min"] = best_known[instance_id_lookup]
    
    # Now, add the checker computed total objective value.
    result["checker_total_objective"] = extract_checker_value(file_path)

    return result

def main():
    if len(sys.argv) < 2:
        print("Usage: python extract.py <directory_with_log_files>")
        sys.exit(1)
        
    directory = sys.argv[1]
    output_csv = "results.csv"
    
    # Collect log files (adjust file extensions if needed)
    log_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.startswith("X"):
                continue
            if file.endswith(".txt") or file.endswith(".json"):
                log_files.append(os.path.join(root, file))
    
    all_data = []
    for log_file in log_files:
        try:
            data = extract_values(log_file)
            all_data.append(data)
        except Exception as e:
            print(f"Error processing file {log_file}: {e}")
    
    # CSV header including the new checker_total_objective field.
    fields = [
        "log_instance", "allowed_time", "time_to_build_model", "num_variables", 
        "num_constraints", "total_solve_time", "time_to_first_solution", 
        "first_solution_value", "gap_to_optimum", "best_solution_value",
        "best_obj_15min", "checker_total_objective"
    ]
    
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fields)
        writer.writeheader()
        for row in all_data:
            writer.writerow(row)
    
    print(f"CSV file saved as {output_csv}")

if __name__ == "__main__":
    main()
