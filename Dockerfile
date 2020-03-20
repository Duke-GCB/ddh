FROM rocker/tidyverse:3.6.1
RUN install2.r -r https://cloud.r-project.org --deps TRUE here janitor corrr beepr enrichR moderndive pander vroom rentrez feather optparse tidytext widyr
