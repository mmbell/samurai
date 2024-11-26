import argparse
import netCDF4 as nc

def main(input_file, output_file, variables_to_keep):

    # Open the original NetCDF file
    with nc.Dataset(input_file, 'r') as src, nc.Dataset(output_file, 'w') as dst:
        # Copy global attributes
        dst.setncatts({attr: src.getncattr(attr) for attr in src.ncattrs()})
        
        # Copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
        
        # Copy selected variables
        for name, variable in src.variables.items():
            if name in variables_to_keep:
                new_var = dst.createVariable(name, variable.datatype, variable.dimensions)
                new_var.setncatts({attr: variable.getncattr(attr) for attr in variable.ncattrs()})
                new_var[:] = variable[:]

    print(f"New NetCDF file created at {output_file} with only the selected variables.")

if __name__ == "__main__":
    # Define the default variables to keep
    default_variables = ['U', 'V', 'W', 'DIV', 'VORT']

    # Create a string representation of the default variables
    default_variables_str = ', '.join(default_variables)

    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Reduce NetCDF file by keeping only specified variables.')
    parser.add_argument('input_file', type=str, help='Path to the input NetCDF file.')
    parser.add_argument('output_file', type=str, help='Path to the output NetCDF file.')
    parser.add_argument('--variables', type=str, nargs='+', default=default_variables,
                        help=f'List of variables to keep (default: {default_variables_str}).')

    # Parse the arguments
    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args.input_file, args.output_file, args.variables)