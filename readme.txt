This repository contains the required R code and data files to produce the Mo-U mass balance model plots (Figures 2 and 3a) presented in ‘Persistent global marine euxinia in the early Silurian’ 
Richard G. Stockey, Devon B. Cole, Noah J. Planavsky, David K. Loydell, Jiří Frýda and Erik A. Sperling

Before running the R scripts included, download this folder and set it as your R working directory. Required data files should then load when called in each script, and pdf files will save within the same folder.

To replicate the analyses and plots presented here, the following R packages are required:
deeptime
deSolve
dplyr
egg
fANCOVA
FME
ggplot2
ggridges
MASS
plyr
RColorBrewer
readr
reshape2

The deeptime package is currently only available on GitHub, and will need to be installed from there to reproduce Fig. 3a. This can be achieved by removing # symbols and running the R commands on lines 15-16 of Fig. 3a script.

All other packages can be installed from CRAN. These scripts have been tested using R version 3.5.0 - Copyright (C) 2018 The R Foundation for Statistical Computing.
