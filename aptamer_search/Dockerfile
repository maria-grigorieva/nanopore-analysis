FROM continuumio/miniconda3:4.12.0
WORKDIR /app
COPY environment.yml .
RUN conda env create --file environment.yml
ENV PATH /opt/conda/envs/nanopore/bin:$PATH
COPY aptamer_search.py .
RUN chmod +x aptamer_search.py
ENTRYPOINT [ "/opt/conda/envs/nanopore/bin/python", "aptamer_search.py" ]