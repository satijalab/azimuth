FROM satijalab/seurat:3.2.0

RUN apt-get update
RUN apt-get install -y libv8-dev

RUN R -e "install.packages('remotes')"

COPY Rprofile.site /usr/lib/R/etc/
COPY . /root/seurat-mapper

RUN R -e "remotes::install_local('/root/seurat-mapper')"

EXPOSE 3838

CMD ["R", "-e", "SeuratMapper::AzimuthApp(reference='/reference-data')"]
