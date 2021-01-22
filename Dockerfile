FROM satijalab/seurat:3.2.3

RUN apt-get update
RUN apt-get install -y libv8-dev

RUN mkdir lzf 
WORKDIR /lzf
RUN wget https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf_filter.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf_filter.h 
RUN mkdir lzf 
WORKDIR /lzf/lzf
RUN wget https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf_c.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf_d.c https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzfP.h https://raw.githubusercontent.com/h5py/h5py/3.0.0/lzf/lzf/lzf.h
WORKDIR /lzf
RUN gcc -O2 -fPIC -shared lzf/*.c lzf_filter.c $(pkg-config --cflags --libs hdf5) -o liblzf_filter.so
WORKDIR /
ENV HDF5_PLUGIN_PATH=/lzf

RUN R --no-echo -e "install.packages('remotes')"

COPY Rprofile.site /usr/lib/R/etc/
RUN R --no-echo -e "install.packages(c('DT', 'future', 'ggplot2',  'googlesheets4', 'hdf5r', 'htmltools', 'httr', 'patchwork', 'rlang', 'shiny', 'shinyBS', 'shinydashboard', 'shinyjs', 'stringr', 'withr', 'BiocManager'), repo='https://cloud.r-project.org')"
RUN R --no-echo -e "remotes::install_github(c('immunogenomics/presto', 'jlmelville/uwot', 'mojaveazure/seurat-disk', 'satijalab/seurat@feat/descartes'))"
RUN R --no-echo -e "BiocManager::install('glmGamPoi')"

ARG SEURAT_VER=unknown
RUN echo "$SEURAT_VER"
RUN R --no-echo -e "remotes::install_github('satijalab/seurat@feat/descartes')"

ARG AZIMUTH_VER=unknown
RUN echo "$AZIMUTH_VER"
COPY . /root/azimuth
RUN R --no-echo -e "remotes::install_deps('/root/azimuth')"
RUN R --no-echo -e "install.packages('/root/azimuth', repos = NULL, type = 'source')"

EXPOSE 3838

CMD ["R", "-e", "Azimuth::AzimuthApp(reference='/reference-data')"]
