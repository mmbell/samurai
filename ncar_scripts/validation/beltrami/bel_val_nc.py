import numpy as np
import xarray as xr
import argparse
import os
import sys

# ANSI escape codes for colors (Blue is good, Orange is bad)
BLUE = '\033[94m'
ORANGE = '\033[93m'
RESET = '\033[0m'


# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Usage: python3 bel_val_nc.py ref_analysis.nc test_analysis.nc')
parser.add_argument('ref', type=str, help='Reference analytical solution netcdf file')
parser.add_argument('run', type=str, help='Test analysis netcdf file')

args = parser.parse_args()

ref = args.ref
run = args.run

if ref is None or run is None:
    raise ValueError('Please provide reference and run log files. \n \
                      Usage: python3 bel_val_nc.py ref_analysis.nc test_analysis.nc ')

# Pull specified variables from the files
variables = ['U', 'V', 'W', 'VORT', 'DIV']

# Error codes
# 0: No error
# 1: File not found
# 2: Variable not found
# 3: RMSE exceeds threshold
error = 0

error_messages = {
    0: "Validation successful. All RMSE values are within acceptable thresholds.",
    1: "Error 1: File not found",
    2: "Error 2: RMSE exceeds threshold"
}

# Define RMSE thresholds for each variable
rmse_thresholds = {
    'U': 0.03396,
    'V': 0.03396,
    'W': 0.03957,
    'VORT': 7.132e-5,
    'DIV': 8.618e-5
}

# Opens the netcdf file and returns only the variables specified
def extract_variables(file_name, variables):
    # Read the file into an xarray Dataset
    ds = xr.open_dataset(file_name)

    # Select the specified variables
    ds = ds[variables]

    return ds

# Calculate the root mean square error for the given variable 
# between two datasets
def calculate_rmse(data1, data2, variable):
    # Extract the variable from both datasets
    var1 = data1[variable].values
    var2 = data2[variable].values

    # Calculate the square of the difference between the two variables
    diff_sq = (var1 - var2) ** 2

    # Take the mean of the squared differences
    mean_diff_sq = np.mean(diff_sq)

    # Take the square root of the mean to get the RMSE
    rmse = np.sqrt(mean_diff_sq)

    return rmse


if os.path.isfile(ref):
    print("Reference file exists: ", ref)
else:
    print("Reference file does not exist at: ", ref)
    error=1

if os.path.isfile(run):
    print("Comparison file exists: ", run)
else:
    print("Comparison file does not exist at: ", run)
    error=1

#exit on file not found error
if error==1:
    print(error_messages[error])
    sys.exit(error)

print("Extracting variables: ", variables)
refVars = extract_variables(ref, variables)
runVars = extract_variables(run, variables)

for variable in variables:
    # Calculate the RMSE for this variable
    rmse = calculate_rmse(refVars, runVars, variable)

    # Check if the RMSE exceeds the threshold
    if rmse > rmse_thresholds[variable]:
        print(f'{ORANGE}RMSE for {variable} exceeds the threshold of {rmse_thresholds[variable]}: {rmse:.5f}{RESET}')
        error = 2
    else:
        print(f'{BLUE}RMSE for {variable}: {rmse:.5f}{RESET}')

if error in error_messages:
    print(error_messages[error])
else:
    print("Unknown error")

sys.exit(error)
