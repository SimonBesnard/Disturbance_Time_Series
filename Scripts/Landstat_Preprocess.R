# Script to pre-process Landsat TS analysis
## Author: Simon Besnard
## 01.06.2015
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
rasterOptions(progress = 'text') # to show progress bar in raster calculations

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

#Number of observations
obs <- countObs(s)
plot(obs)
summary(obs)

# valid observations
obs <- countObs(s, as.perc=TRUE)
summary(s)

# % NA per pixel
percNA <- 100 - countObs(s, as.perc=TRUE)
plot(percNA, main="percent NA per pixel")

#Compute mean VI
meanVI <- summaryBrick(s, fun=mean, na.rm=T)
plot(meanVI)

# define a function that takes a vector as an argument
checkThresh <- function(x){
  # first, get rid of NA's
  x <- x[!is.na(x)]
  # if there still values left, count how many are above the threshold
  # otherwise, return a 0
  if(length(x) > 0){
    y <- length(x[x > 7000])
  } else {
    y <- 0
  }
  # return the value
  return(y)
}

# pass this functiom to summaryBrick
customStat <- summaryBrick(s, fun=checkThresh)
plot(customStat, main = "# of observations where NDVI > 0.7")

# median values for all layers
medVI <- summaryBrick(s, fun=median, na.rm=TRUE)
# only ETM+ layers
medVI_ETM <- summaryBrick(s, fun=median, na.rm=TRUE, sensor="ETM+")
# all layers between 1984 and 2005 (inclusive)
medVI_84_00 <- summaryBrick(s, fun=median, na.rm=TRUE, minDate="1984-06-21", 
                            maxDate="2000-12-31")
# all layers after 2005
medVI_01_14 <- summaryBrick(s, fun=median, na.rm=TRUE, minDate=c(2001, 1))

# plot and compare
op <- par(mfrow=c(2, 2))
plot(medVI, main = "median NDMI")
plot(medVI_ETM, main = "only ETM+")
plot(medVI_84_00, main = "1984-2000")
plot(medVI_01_14, main = "2001-2014")

# Compute annual mean 
annualMed <- annualSummary(s, fun=median, na.rm=TRUE)
plot(annualMed)

# 2.Use of Spatial BFASTMonitor

# run bfmPixel() in interactive mode with a monitoring period 
plot(s,50)
bfm <- bfmPixel(s, start=c(2000, 1), interactive=TRUE)
plot(bfm$bfm)

# run bfmPixel() with target cell 
targcell <- 3592
bfm <- bfmPixel(s, cell=targcell, start=c(2000, 1))

# use a harmonic model only
bfm1 <- bfmPixel(s, cell=targcell, start=c(2000, 1), 
                 formula=response~harmon, plot=TRUE)

# same, but with an order of 1
bfm2 <- bfmPixel(s, cell=targcell, start=c(2000, 1), 
                 formula=response~harmon, order=1, plot=TRUE)

# only trend
bfm3 <- bfmPixel(s, cell=targcell, start=c(2000, 1), 
                 formula=response~trend, plot=TRUE)

# bfmPixel using a 2-years monitoring period
bfm4 <- bfmPixel(s, cell=targcell, start=c(1999, 1), 
                 monend=c(2005, 1), plot=TRUE)

# apply bfmPixel only on ETM+ data
bfm5 <- bfmPixel(s, cell=targcell, start=c(2000, 1), 
                 sensor="ETM+", plot=TRUE)

# 3. Use of bfmSpatial over flux sites

# run bfmSpatial over flux sites for a monitoring period of 1999 with multicore
t1 <- system.time(
  bfm_Au_Tum <- bfmSpatial(s, start=c(1999, 1), order=1, monend=c(2005, 1), mc.cores=2)
)

# Save bfm outuput in RDS file
saveRDS(bfm_Au_Tum, 'Output/bfm_Au_Tum.rds')

# extract change raster
change <- raster(bfm, 1)
months <- changeMonth(change)
# set up labels and colourmap for months
monthlabs <- c("jan", "feb", "mar", "apr", "may", "jun", 
               "jul", "aug", "sep", "oct", "nov", "dec")
cols <- rainbow(12)
plot(months, col=cols, breaks=c(1:12), legend=FALSE)
# insert custom legend
legend("bottomright", legend=monthlabs, cex=0.5, fill=cols, ncol=2)

# extract magn raster
magn <- raster(bfm, 2)
# make a version showing only breakpoing pixels
magn_bkp <- magn
magn_bkp[is.na(change)] <- NA
op <- par(mfrow=c(1, 2))
plot(magn_bkp, main="Magnitude: breakpoints")
plot(magn, main="Magnitude: all pixels")

# Run bfmSpatial using a 1-year monitoring period
par(op)
magn <- raster(bfm, 2) / 10000
bfm99 <- bfmSpatial(s, start=c(1999, 1), monend=c(2000, 1), order=1)

# extract change
change99 <- raster(bfm99, 1)
# extract and rescale magnitude
magn99 <- raster(bfm99, 2) / 10000
# remove all non-breakpoints
magn99[is.na(change99)] <- NA

# extract and rescale magnitude and apply a -500 threshold
magn99thresh <- magn99
magn99thresh[magn99 > -0.05] <- NA
# compare
op <- par(mfrow=c(1, 2))
plot(magn99, main="magnitude")
plot(magn99thresh, main="magnitude < -0.05")

par(op)

magn99_sieve <- areaSieve(magn99thresh, thresh=1800)
magn99_areasieve <- areaSieve(magn99thresh)
magn99_as_rook <- areaSieve(magn99thresh, directions=4)

# compare all magn rasters
op <- par(mfrow=c(2, 2))
plot(magn99thresh, main="magnitude")
plot(magn99_sieve, main="pixel sieve")
plot(magn99_areasieve, main="0.5ha sieve")
plot(magn99_as_rook, main="0.5ha sieve, rook's case")

par(op)

changeSize_queen <- clumpSize(magn99_areasieve)
changeSize_rook <- clumpSize(magn99_areasieve, directions=4)
# compare
op <- par(mfrow=c(1, 2))
plot(changeSize_queen, col=bpy.colors(50), main="Clump size: Queen's case")
plot(changeSize_rook, col=bpy.colors(50), main="Clump size: Rook's case")

par(op)

changeSize <- clumpSize(magn99_areasieve, f=900/10000)
plot(changeSize, col=bpy.colors(50), main="Clump size (hectares)")

changeSize <- clumpSize(magn99_areasieve, f=900/10000, stats=TRUE)
print(changeSize$stats)

# 3. Make plot of the time series with bfastApp

# Generate SpatialPoints object and visualize
sp <- sampleRegular(s, size = 160, sp=TRUE)

# Extract samples and prepare rds file
zooExtract(x = s, sample = sp, file = 'Output/BR_Sa3_Zoo.rds')
runGitHub(repo = 'bfastApp', username = 'dutri001')

# 4. Plot time series using a sample grid

# Source functions
source('Function/bpPhenoShiny.R')
source('Function/ggplot_bfastIR.R')

# Load ts objects (zooExtract output)
ts <- readRDS('Output/BR_SA2_Zoo.rds')
ts <- ts / 10000

# Run functions
bp <- bpPhenoShiny(x = ts[,132], formula = response ~ trend, h = 0.40)
gg <- ggplot.bfastIR(bp, formula = response ~ trend, order = 2, numbering = TRUE)

# Then you can continue customizing the gg object
gg <- gg +
  xlab('Time') +
  ylab('NDMI (-)') +
theme(panel.background = element_rect(fill = 'white'),
      panel.border = element_rect(colour = 'black', fill = NA),
      title = element_text(colour = 'black'),
      axis.text = element_text(colour = 'black'))

gg
ggsave(gg, filename = 'Latex/Figures/BR_Sa2_Pheno.png', width = 9, height = 4)

# 5. Plot BfastSpatial outputs

#Create dataframe for ggplot
ggr<- fortify(bfm_BR_Sa2)
ggr[,1][ggr[,1] == -Inf] <- NA
ggr[,1]<-date_decimal(ggr[,1])
ggr[,1]<- sapply(ggr[,1],  function(x) year(x))

sp.df <- data.frame(coordinates(sp))

#Create ggplot files
g1<-ggplot() +
  geom_raster(data = ggr, aes(x,y,fill = factor(values.layer.1))) +
  geom_point(data=sp.df, aes(x=x, y=y), color="black", shape= 3, size=1)+
  scale_fill_manual(name= "Year of break detected", 
                      values = c("#E69F00", "#56B4E9", 
                                 "#009E73", "#F0E442", 
                                 "#0072B2","#D55E00", "#66FF22"))+
  coord_equal() +
  scale_x_continuous(name = 'Long', expand=c(0, 0)) +
  scale_y_continuous(name = 'Lat', expand=c(0, 0))+
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = 'black', fill = NA),
        title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'))

g2<-ggplot() +
  geom_raster(data = ggr, aes(x,y,fill = values.layer.2/10000)) +
  geom_point(data=sp.df, aes(x=x, y=y), color="black", shape= 3, size=1)+
  scale_fill_gradient2(low = "red", mid = "white",
                        high = "green", midpoint = 0, space = "Lab",
                        na.value = "grey50", guide = "colourbar", name= "Magnitude of disturbance", limits=c(-1,1))+
  coord_equal() +
  scale_x_continuous(name = 'Long', expand=c(0, 0)) +
  scale_y_continuous(name = 'Lat', expand=c(0, 0))+
  theme(panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = 'black', fill = NA),
        title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'))

# Plot two graphs together
gA <- ggplotGrob(g1)
gB <- ggplotGrob(g2)
maxWidth = grid::unit.pmax(gA$widths[1:2], gB$widths[1:2])
gA$widths[1:2] <- as.list(maxWidth)
gB$widths[1:2] <- as.list(maxWidth)
g3 <- arrangeGrob(
  gA, gB, nrow = 1, heights = c(2, 2))

# Save plots
png(filename="Latex/Figures/BR_Sa2_Bfm.png", width=12, height=3, units="in", res=300)
g3
dev.off()