FROM rocker/r-ver:3.5

RUN apt update && apt install -y --no-install-recommends kmc bedtools ncbi-blast+ libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev pandoc gawk

RUN install2.r tidyverse lsa DT plotly rmarkdown optparse

WORKDIR /app
