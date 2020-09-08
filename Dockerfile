FROM satijalab/seurat:3.2.0

RUN apt-get update
RUN apt-get install -y libv8-dev

RUN mkdir lzf 
WORKDIR /lzf
RUN wget https://raw.githubusercontent.com/h5py/h5py/master/lzf/lzf_filter.c https://raw.githubusercontent.com/h5py/h5py/master/lzf/lzf_filter.h
RUN mkdir lzf 
WORKDIR /lzf/lzf
RUN wget https://raw.githubusercontent.com/h5py/h5py/master/lzf/lzf/lzf_c.c https://raw.githubusercontent.com/h5py/h5py/master/lzf/lzf/lzf_d.c https://raw.githubusercontent.com/h5py/h5py/master/lzf/lzf/lzf.h https://raw.githubusercontent.com/h5py/h5py/master/lzf/lzf/lzfP.h
WORKDIR /lzf
RUN gcc -O2 -fPIC -shared lzf/*.c lzf_filter.c $(pkg-config --cflags --libs hdf5) -o liblzf_filter.so
WORKDIR /
ENV HDF5_PLUGIN_PATH=/lzf

RUN R -e "install.packages('remotes')"

RUN wget https://raw.githubusercontent.com/pshved/timeout/master/timeout && chmod +x timeout

COPY Rprofile.site /usr/lib/R/etc/
COPY . /root/seurat-mapper

RUN R -e "remotes::install_local('/root/seurat-mapper')"

EXPOSE 3838

CMD ["R", "-e", "SeuratMapper::AzimuthApp(reference='/reference-data')"]
