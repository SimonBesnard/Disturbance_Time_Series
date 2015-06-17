# Script to pre-process MODIS TS analysis
## Author: Simon Besnard
## 10.06.2015
###################################
## Load the necessary packages
# install.packages('devtools')
# install_github('dutri001/bfastSpatial')
# install.packages("rgdal")
# install.packages("raster")
# install.packages("gdalUtils")
# install.packages("plyr")
# install.packages("rgeos")
# install.packages("reshape")
# install.packages('maptools')
# install.packages('caret')
# install.packages("igraph")
# install.packages('shiny')
# install.packages('lubridate')
# devtools::install_github('dutri001/ggSpatial')
# install.packages("gridExtra")
# install.packages("xts", repos="http://r-forge.r-project.org")

# load the package
library (sp)
library(devtools)
library(rgdal)
library(bfastSpatial)
library(bfast)
library (raster)
library(rgdal)
library(gdalUtils)
library (rgeos)
library(plyr)
library (caret)
library(parallel)
library (stringr)
library(shiny)
library(lubridate)
library(ggplot2)
library (strucchange)
library (ggSpatial)
library(gridExtra)
library (RNetCDF)
library(xts)
rasterOptions(progress = 'text') # to show progress bar in raster calculations

# 1.Load MODIS data 
MODIS_AU_Tum<-open.nc (file.path(path, 'Modis_Data/BR-Sa3.MODIS.TERRA.C5.VIs.JCP1.2001-2014.daily.nc'), write=T)
df<-read.nc(MODIS_AU_Tum)

# 2. Create a time series
Au_Tum_Zoo<- na.omit(zoo(df$NDVI_250_original, order.by = as.Date(df$time, origin=as.Date("1582-10-15"))))
Au_Tum_Zoo <- zoo(Au_Tum_Zoo, index(Au_Tum_Zoo))

# 3.Plot time series and detect breaks

# Source functions
source('Function/bpPhenoShiny.R')
source('Function/ggplot_bfastIR.R')
bp <- bpPhenoShiny(x = Au_Tum_Zoo, formula = response ~ trend + harmon, order=2,  breaks=1, h=0.25)
gg <- ggplot.bfastIR(bp, formula = response ~ trend + harmon, order = 2, numbering = TRUE)

# Then you can continue customizing the gg object
gg <- gg +
  xlab('Time') +
  ylab('NDVI (-)') +
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = 'black', fill = NA),
        title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'))

