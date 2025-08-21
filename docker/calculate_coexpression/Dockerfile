# syntax=docker/dockerfile:1

FROM rocker/r-ver:4.4.1 

RUN R -q -e "install.packages(c('remotes', 'data.table'))"
RUN R -q -e "remotes::install_github('PavlidisLab/aggtools')"

# Get the coexpression script and make it executable
COPY scripts/calculate_coexpression.R /usr/local/bin/
RUN chmod +x /usr/local/bin/calculate_coexpression.R
