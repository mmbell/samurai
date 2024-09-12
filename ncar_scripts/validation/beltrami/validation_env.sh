#!/bin/bash

# Set the working directory
cd /app

# Ensure pip is installed
apt-get update && apt-get install --reinstall zlib1g
apt-get install -y python3-venv
python3 -m venv /app/venv && \
    /app/venv/bin/pip install --upgrade pip && \
    /app/venv/bin/pip install numpy xarray netCDF4 scipy cftime 

# Make the validation script executable
chmod +x /app/samurai/ncar_scripts/validation/beltrami/bel_val_nc.py

echo "#!/bin/bash" > /app/run_validation.sh
echo "cd /app/samurai/build/release/bin/" >> /app/run_validation.sh
echo "./samurai -params /app/samurai/ncar_scripts/TDRP/beltrami.tdrp | tee /app/samurai/ncar_scripts/validation/beltrami/test_log" >> /app/run_validation.sh
echo "cd /app/samurai/ncar_scripts/validation/beltrami/" >> /app/run_validation.sh
echo "/app/venv/bin/python3 bel_val_nc.py samurai_XYZ_wind_analysis_ref.nc /app/samurai/build/release/bin/samurai_XYZ_wind_analysis.nc" >> /app/run_validation.sh

chmod +x /app/run_validation.sh

# Run the validation script
#./app/run_validation.sh
