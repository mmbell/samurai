import re
import argparse
import math
import numpy as np

# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Usage: python validation.py refernce_log_file run_log_file')
parser.add_argument('ref', type=str, help='Reference log file')
parser.add_argument('run', type=str, help='Run log file')

args = parser.parse_args()

ref = args.ref
run = args.run

if ref is None or run is None:
    raise ValueError('Please provide reference and run log files')

def parseFile(file):
    pattern = r"CG iteration\s+(\d+):\s+r_norm =\s+(\d+(\.\d+)?)\s+rel_resid =\s+(\d+(\.\d+)?)"

    with open(file) as f:
        data = f.read()

    matches = re.findall(pattern, data)

    # Initialize lists
    cg_iterations = []
    r_norms = []
    rel_resids = []

    # Populate lists
    for match in matches:
        cg_iterations.append(int(match[0]))
        r_norms.append(float(match[1]))
        rel_resids.append(float(match[3]))

    # Create dictionary
    data_dict = {
        'CG iteration': cg_iterations,
        'r_norm': r_norms,
        'rel_resid': rel_resids
    }

    return data_dict

def calculate_rms(dict1, dict2):
    if set(dict1.keys()) != set(dict2.keys()):
        raise ValueError("Dictionaries have different keys")

    rms_values = {}
    for key in dict1.keys():
        diff = np.array(dict1[key]) - np.array(dict2[key])
        square_diff = np.square(diff)
        mean_square_diff = np.mean(square_diff)
        rms = np.sqrt(mean_square_diff)
        rms_values[key] = rms

    return rms_values

ref_data = parseFile(ref)
run_data = parseFile(run)
rms_values = calculate_rms(ref_data, run_data)
for key, rms in rms_values.items():
    print(f'RMS value for {key}: {rms}')