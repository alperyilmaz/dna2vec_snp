FROM rocker/r-ver:3.5

RUN apt update && apt install -y --no-install-recommends jellyfish bedtools ncbi-blast+ libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev

RUN install2.r tidyverse lsa DT plotly

