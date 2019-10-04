FROM rocker/tidyverse:3.6.1
RUN install2.r --deps TRUE here janitor corrr beepr enrichR moderndive pander vroom
RUN apt-get install texlive-xetex -y
ADD ./code /depmap/code
ADD ./.Renviron /depmap/.Renviron
RUN mkdir /depmap/data
# render method creates temp files in /depmap/code
RUN chmod a+w /depmap/code
WORKDIR /depmap
CMD Rscript /depmap/code/correlate_data.R
