# Author: Supreeth Suresh
# Affiliation: National Center for Atmospheric Research, Computational and Information Systems Laboratory
# Created: 03-29-2024

FROM ssuresh1809/samurai:nvhpc24_3

LABEL maintainer="Supreeth Suresh <ssuresh@ucar.edu>"

RUN mkdir -p /app/

WORKDIR /app/

COPY . /app/samurai

RUN cd samurai && \
    export SAMURAI_ROOT=$(pwd) && \
    export LROSE_INSTALL_DIR=/usr/local && \
    mkdir build && \
    cd build && \
    cmake -DUSE_GPU=false -DDEBUG_COMPILE=false .. && \
    make -j 2 VERBOSE=1

RUN sed -i 's|data_directory = "/glade/campaign/cisl/asap/samurai/data/beltrami";|data_directory = "/app/data/samurai_data/beltrami";|g' /app/samurai/ncar_scripts/TDRP/beltrami.tdrp

ENV PATH="/app/samurai/build:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
ENV PATH="/usr/local/bin:${PATH}"
#RUN pip install --no-cache-dir numpy

RUN echo "#!/bin/bash \n\
cd /app/samurai/build/release/bin/ \n\
./samurai -params /app/samurai/ncar_scripts/TDRP/beltrami.tdrp | tee /app/samurai/ncar_scripts/validation/beltrami/test_log \n\
" > /app/run_samurai.sh
WORKDIR /app/
RUN chmod +x /app/run_samurai.sh
#RUN chmod +x /app/samurai/ncar_scripts/validation/beltrami/validation_env.sh
#RUN chmod +x /app/samurai/ncar_scripts/validation/beltrami/bel_val_nc.py

CMD [ "./run_samurai.sh" ]

#CMD ["python3", "/app/samurai/ncar_scripts/validation/beltrami/bel_val_nc.py", "/app/samurai/ncar_scripts/validation/beltrami/samurai_XYZ_wind_analysis_ref.nc", "/app/samurai/build/release/bin/samurai_XYZ_wind_analysis.nc"]
#python3 /app/samurai/ncar_scripts/validation/beltrami/bel_val_nc.py /app/samurai/ncar_scripts/validation/beltrami/samurai_XYZ_wind_analysis_ref.nc /app/samurai/build/release/bin/samurai_XYZ_wind_analysis.nc > /app/beltrami_output.log \n\
#python /app/samurai/ncar_scripts/validation/beltrami/compare_log_file.py /app/samurai/ncar_scripts/validation/beltrami/ref_log /app/samurai/ncar_scripts/validation/beltrami/test_log  \n\
