FROM satijalab/azimuth:latest

RUN R --no-echo -e "BiocManager::install('cowplot')"


ARG AZIMUTH_VER=unknown
RUN echo "$AZIMUTH_VER"
COPY . /root/azimuth
RUN R --no-echo -e "install.packages('/root/azimuth', repos = NULL, type = 'source')"

EXPOSE 3838

CMD ["R", "-e", "Azimuth::AzimuthApp(reference='/reference-data')"]
