#!/bin/bash

# Set the working directory
cd /app

# Copy necessary output for validation and delete rest to make space
cp /app/samurai/build/release/bin/samurai_XYZ_wind_analysis.nc /app/samurai/ncar_scripts/validation/beltrami/samurai_XYZ_wind_analysis_ref.nc
rm -rf /opt/nvidia
rm -rf /app/samurai/build
rm -rf /app/samurai/src
rm -rf /app/data

# Ensure pip is installed
apt-get update && \
apt-get install --reinstall -y --no-install-recommends zlib1g python3-venv && \
apt-get clean && \
rm -rf /var/lib/apt/lists/*

# Create a virtual environment and install necessary packages
python3 -m venv /app/venv 

# Upgrade pip
/app/venv/bin/pip install --upgrade pip

# Install packages individually to reduce peak memory usage
/app/venv/bin/pip install numpy 
/app/venv/bin/pip install xarray 
/app/venv/bin/pip install netCDF4 
/app/venv/bin/pip install scipy 
/app/venv/bin/pip install cftime

# Clean up pip cache to free up space
rm -rf /root/.cache/pip

# Make the validation script executable
chmod +x /app/samurai/ncar_scripts/validation/beltrami/bel_val_nc.py

echo "#!/bin/bash" > /app/run_validation.sh
echo "cd /app/samurai/ncar_scripts/validation/beltrami/" >> /app/run_validation.sh
echo "/app/venv/bin/python3 bel_val_nc.py samurai_XYZ_wind_analysis_ref.nc samurai_XYZ_wind_analysis.nc" >> /app/run_validation.sh

chmod +x /app/run_validation.sh

# Run the validation script
#./app/run_validation.sh
