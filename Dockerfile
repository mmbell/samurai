# Author: Supreeth Suresh
# Affiliation: National Center for Atmospheric Research, Computational and Information Systems Laboratory
# Created: 03-29-2024

FROM ssuresh1809/samurai:latest

LABEL maintainer="Supreeth Suresh <ssuresh@ucar.edu>"

RUN mkdir -p /app/
WORKDIR /app/
RUN git clone https://github.com/mmbell/samurai.git && \
    cd samurai && \
    sed -i 's/MODE GPU/MODE CPU/g' CMakeLists.txt && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 2

ENV PATH="/app/samurai/build:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH}"
ENV PATH="/usr/local/bin:${PATH}"

RUN echo "#!/bin/bash \n\
cd /app/samurai/build/release/bin/ \n\
./samurai /app/data/beltrami_preprocessed.xml\
" > /app/run_samurai.sh
WORKDIR /app/
RUN chmod +x /app/run_samurai.sh

CMD [ "./run_samurai.sh" ]