# Script for greenbrown - land surface phenology and trend analysis
## Author: Simon Besnard
## 09.06.2015
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
# install.packages("Kendall")

# load the package
library (sp)
library(rgdal)
library(bfast)
library (raster)
library(rgdal)
library(plyr)
library(lubridate)
library(ggplot2)
library (strucchange)
library(Kendall)
rasterOptions(progress = 'text') # to show progress bar in raster calculations

# Load the greenbrown package
setwd("/home/simonbesnard/R/x86_64-pc-linux-gnu-library/3.2/greenbrown/R/")
files <- list.files(pattern=".R")
for (i in 1:length(files)) source(files[i])

# 1. Pre-process Landsat time series

# Set the location of output and intermediary directories (everything in tmpdir in that case)
srdir <- dirout <- file.path(dirname(rasterTmpFile()), 'bfmspatial') # We use dirname(rasterTmpFile()) instead of rasterOptions()$tmpdir to reduce verbose
dir.create(dirout, showWarning=FALSE)

# Get the directory where the Landsat archives are stored
dir <- file.path(path, 'Landsat_Data/BR_Sa3')
list <- list.files(dir, pattern=glob2rx('*.tar.gz'), full.names=TRUE)

# Run the batch line
processLandsatBatch(x=dir, pattern=glob2rx('*.tar.gz'), outdir=dirout, srdir=srdir, delete=TRUE, vi='ndmi', mask='fmask', keep=0, overwrite=TRUE)

# Visualize one of the layers produced
list <- list.files(dirout, pattern=glob2rx('*.grd'), full.names=TRUE)

# Create a new subdirectory in the temporary directory
dirout <- file.path(dirname(rasterTmpFile()), 'stack')
dir.create(dirout, showWarnings=FALSE)

# Generate a file name for the output stack
stackName <- file.path(dirout, 'stack_BR_Sa2.grd')

# Stack the layers
s <- timeStack(x=list, filename=stackName, datatype='INT2S', overwrite=TRUE)

#2. Analyze trends and trend changes on time series
#Create a time series dataframe
BR_SA2_Zoo<- readRDS('Au_Tum.rds')
BR_SA2_Zoo<- BR_SA2_Zoo/10000
dates<- time(BR_SA2_Zoo)
test <- bfastts(as.vector(BR_SA2_Zoo[,77]), dates, type = 'irregular')

# calculate trend (default method: TrendAAT)
trd <- Trend(test)
plot(trd)

# calculate trend but consider breakpoints
trd <- Trend(test, mosum.pval=0.5, breaks=1, sample.min.length =0.5 )
plot(trd) 

#3. Analyse seasonality
Seasonality(test, plot=TRUE)

#4. Analyze phenology 
spl.trs <- Phenology(test, tsgf="TSGFspline", method="Beck", approach="White",check.seasonality = 2)
plot(spl.trs)
plot(spl.trs$series)

#6.Filling of permanent gaps in time series: FillPermanentGap

# set NA values into winter months to simulate gaps
winter <- (1:length(test))[cycle(test) == 1 | cycle(test) == 2 | cycle(test) == 12]
gaps <- sample(winter, length(winter)*0.3)
ndvi2 <- test
ndvi2[gaps] <- NA
plot(ndvi2)

# check and fill permanent gaps
IsPermanentGap(ndvi2)
fill <- FillPermanentGaps(ndvi2) # default fills with the minimum value
plot(fill, col="red"); lines(test)

#7. Temporal smoothing and gap filling: TsPP

# introduce random gaps 
gaps <- test
gaps[runif(100, 1, length(test))] <- NA

# smoothing and gap filling
tsgf <- TSGFspline(gaps)
plot(gaps, main="Temporal smoothing and gap-filling using splines")
lines(tsgf, col="red")

#8. Detection of phenology events: PhenoDeriv and PhenoTrs

# time series pre-processing
x <- TsPP(tsgf, interpolate=F)[1:1096]

# calculate phenology metrics for first year
PhenoTrs(x, plot=TRUE, approach = "White")

  