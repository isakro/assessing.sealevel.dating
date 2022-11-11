# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/verse:4.2.1

# required
MAINTAINER Isak Roalkvam <isak.roalkvam@iakh.uio.no>

COPY . /assessing.sealevel.dating

RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")'
RUN Rscript -e 'BiocManager::install("IRanges")'

# go into the repo directory
RUN . /etc/environment \
  # Install linux depedendencies here
  # e.g. need this for ggforce::geom_sina
  && sudo apt-get update \
  && sudo apt-get install libudunits2-dev -y \
  && apt-get install gdal-bin -y \
  # build this compendium package
  && R -e "devtools::install('/assessing.sealevel.dating', dep=TRUE)" \
  # && R -e "source('/assessing.sealevel.dating/analysis/script/00dtm_prep.R')"
  && R -e "source('/assessing.sealevel.dating/analysis/script/01site_prep.R')"
  # && R -e "source('/assessing.sealevel.dating/analysis/script/02displacement_prep.R')"
  # && R -e "source('/assessing.sealevel.dating/analysis/script/03dtm_edit.R')"
  # && R -e "source('/assessing.sealevel.dating/analysis/script/05site_analysis.R')"
  # && R -e "source('/assessing.sealevel.dating/analysis/script/06simulation_results.R')"
  # && R -e "source('/assessing.sealevel.dating/analysis/script/07shoreline_date.R')"
  # && R -e "source('/assessing.sealevel.dating/analysis/script/08backwards_dating.R')"
  # && R -e "source('/assessing.sealevel.dating/analysis/script/09redate.R')"
  # render the manuscript into a docx, you'll need to edit this if you've
  # customised the location and name of your main Rmd file
  && R -e "rmarkdown::render('/assessing.sealevel.dating/analysis/paper/paper.Rmd')"

