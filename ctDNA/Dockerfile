FROM fanyucai1/bio_tumor:latest
RUN mkdir -p /script/core && mkdir -p /software/CNVkit/cnvkit-master/data/ && \
    mkdir -p /software/genefuse/genes/ && mkdir -p /data/Database/hg19/ && \
    mkdir -p /data/Database/CNV_baseline/
COPY script/ /script/
COPY refFlat.txt /software/CNVkit/cnvkit-master/data/
COPY access-5k-mappable.hg19.bed /software/CNVkit/cnvkit-master/data/
COPY cancer.hg19.csv  /software/genefuse/genes/
COPY druggable.hg19.csv /software/genefuse/genes/
RUN chmod 777 -R /script/
ENV PATH /software/python3/Python-v3.7.0/bin/:$PATH