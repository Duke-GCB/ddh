FROM rocker/tidyverse:3.6.1
RUN install2.r --deps TRUE here janitor corrr beepr enrichR moderndive pander vroom
ADD ./code /depmap/code
ADD ./.Renviron /depmap/.Renviron
RUN mkdir /depmap/data
WORKDIR /depmap
CMD Rscript /depmap/code/correlate_data.R
