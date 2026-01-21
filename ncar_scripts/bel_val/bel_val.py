import numpy as np
import pandas as pd
import os


def extract_variables(file_name, variables):
    # Read the file into a pandas DataFrame
    df = pd.read_csv(file_name, sep="\t", comment="#", header=0, usecols=variables)

    # Convert the DataFrame to an xarray Dataset
    #ds = xr.Dataset.from_dataframe(df)

    return df

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

# +

# Specify the paths of the two output files - one from the reference simulation and one from the comparison simulation
ref_path = "/glade/derecho/scratch/cmille73/APAR/samurai/samurai_XYZ_analysis.out"
comp_path = "/glade/derecho/scratch/cmille73/APAR/samurai/samurai_XYZ_analysis.test.out"


if os.path.isfile(ref_path):
    print("Reference file exists: ", ref_path)
else:
    print("Reference file does not exist at: ", ref_path)

if os.path.isfile(comp_path):
    print("Comparison file exists: ", comp_path)
else:
    print("Comparison file does not exist at: ", comp_path)



# +


#pull specified variables from the files
variables = ['u', 'v', 'w', 'Vorticity', 'Divergence']
print("Extracting variables: ", variables)
ref = extract_variables(ref_path, variables)
comp = extract_variables(comp_path, variables)
# -


for variable in variables:
    # Calculate the RMSE for this variable
    rmse = calculate_rmse(ref, comp, variable)

    # Print the RMSE with a nice label
    print(f'RMSE for {variable}: {rmse:.5f}')

