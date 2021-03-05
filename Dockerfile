FROM satijalab/azimuth:latest

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
