import numpy as np
import xarray as xr
import argparse
import os
import sys


# Set up command line argument parsing
parser = argparse.ArgumentParser(description='Usage: python3 bel_val_nc.py ref_analysis.nc test_analysis.nc')
parser.add_argument('ref', type=str, help='Reference analytical solution netcdf file')
parser.add_argument('run', type=str, help='Test analysis netcdf file')

args = parser.parse_args()

ref = args.ref
run = args.run

if ref is None or run is None:
    raise ValueError('Please provide reference and run log files')

# Pull specified variables from the files
variables = ['u', 'v', 'w', 'Vorticity', 'Divergence']

# Define RMSE thresholds for each variable
rmse_thresholds = {
    'u': 0.03396,
    'v': 0.03396,
    'w': 0.03957,
    'Vorticity': 7.132e-5,
    'Divergence': 8.618e-5
}

def extract_variables(file_name, variables):
    # Read the file into an xarray Dataset
    ds = xr.open_dataset(file_name)

    # Select the specified variables
    ds = ds[variables]

    return ds

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

if os.path.isfile(run):
    print("Comparison file exists: ", run)
else:
    print("Comparison file does not exist at: ", run)


print("Extracting variables: ", variables)
refVars = extract_variables(ref, variables)
runVars = extract_variables(run, variables)

for variable in variables:
    # Calculate the RMSE for this variable
    rmse = calculate_rmse(refVars, runVars, variable)

    # Print the RMSE with a nice label
    print(f'RMSE for {variable}: {rmse:.5f}')

    # Check if the RMSE exceeds the threshold
    if rmse > rmse_thresholds[variable]:
        print(f'RMSE for {variable} exceeds the threshold of {rmse_thresholds[variable]}')
        sys.exit(0)

print("All RMSE values are within acceptable thresholds.")
