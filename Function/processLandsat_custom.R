processLandsat_custom <- function (x, vi = "ndvi", srdir, outdir, untar = TRUE, delete = FALSE, 
                                   mask = NULL, L = 0.5, ...) 
{
  if (untar) {
    ex <- extension(x)
    if (ex == ".gz") {
      tarlist <- untar(x, list = TRUE)
    }
    else if (ex == ".zip") {
      tarlist <- unzip(x, list = TRUE)$Name
    }
    else {
      stop("The archive is neither tar.gz nor .zip; we don't know what to do with that.")
    }
    if (any(grepl(pattern = "^.*\\.hdf$", x = tarlist))) {
      x0 <- grep(pattern = "^.*\\.hdf$", x = tarlist, value = TRUE)
    }
    else if (any(grepl(pattern = "^.*\\.tif$", x = tarlist))) {
      if (any(grepl(pattern = sprintf("^.*%s\\.tif$", vi), 
                    x = tarlist))) {
        x0 <- grep(pattern = sprintf("^.*%s\\.tif$", 
                                     vi), x = tarlist, value = TRUE)
      }
      else {
        
        if (vi == "ndvi") {
          viFormula <- .ndvi()
        }
        else if (vi == "evi") {
          viFormula <- .evi()
        }
        else if (vi == "nbr") {
          viFormula <- .nbr()
        }
        else if (vi == "savi") {
          viFormula <- .savi(L = L)
        }
        else if (vi == "tcwet") {
          viFormula <- .tcwet()
        }
        else if (vi == "tcbright") {
          viFormula <- .tcbright()
        }
        else if (vi == "tcgreen") {
          viFormula <- .tcgreen()
        }
        else {
          stop("Unsupported vi")
        }
        x0 <- grep(pattern = sprintf("^.*(%s)\\.tif$", 
                                     paste(viFormula$ind, collapse = "|")), x = tarlist, 
                   value = TRUE)
      }
    }
    else {
      stop("Did not find any .tif or .hdf files in the archive")
    }
    if (!is.null(mask)) {
      x0 <- c(x0, grep(pattern = sprintf("^.*%s\\.tif$", 
                                         mask), x = tarlist, value = TRUE))
    }
    if (ex == ".gz") {
      untar(x, files = x0, exdir = srdir)
    }
    else if (ex == ".zip") {
      unzip(x, files = x0, exdir = srdir)
    }
    x <- file.path(srdir, x0)
  }
  name <- str_extract(string = basename(x[1]), "(LT4|LT5|LE7|LC8)\\d{13}")
  #   writeRaster(sr2vi_custom(x = x, vi = vi, filename = sprintf("%s/%s.%s.grd", 
  #                                            outdir, vi, name), datatype = "INT2S", mask = mask),filename = sprintf("%s/%s.%s.grd", outdir, vi, name),overwrite=TRUE)
  sr2vi_custom(x = x, vi = vi, filename = sprintf("%s/%s.%s.grd", outdir, vi, name), datatype = "INT2S", mask = mask)
  if (delete) {
    file.remove(x)
  }
}

for(i in 1:length(list)){
  tc_wet.coef <- c(0.0315,  0.2021,  0.3102,  0.1594, -0.6806, -0.6109)
  tc_bright.coef <- c(0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303)
  tc_green.coef <- c(-0.1603, -0.2819, -0.4934,  0.7940, -0.0002, -0.1446) 
  
  tarlist <- untar(list[i], list = TRUE)
  sr_files <- tarlist[substr(tarlist, regexpr("sr", tarlist), nchar(tarlist)-5) %in% "sr_band"] # select surface reflectance files
  untar(list[i], files=sr_files, exdir=srdir) # extract only metadata file to scene specific output folder
  sr_bands <- lapply(X = list.files(path=srdir, pattern=glob2rx('L*'), recursive=T, full.names=T), FUN = raster)
  s <- stack(sr_bands)
  overlay(s[[1]],s[[2]],s[[3]],s[[4]],s[[5]],s[[6]],fun=function(x1,x2,x3,x4,x5,x7) {x1[]*tc_wet.coef[1] + x2[]*tc_wet.coef[2]+ x3[]*tc_wet.coef[3]+ x4[]*tc_wet.coef[4]+ x5[]*tc_wet.coef[5]+ x7[]*tc_wet.coef[6]},
  filename=paste(out.tcwet,"/","tcwet.",substr(basename(list[i]),1,16),sep=""),overwrite=T)
  overlay(s[[1]],s[[2]],s[[3]],s[[4]],s[[5]],s[[6]],fun=function(x1,x2,x3,x4,x5,x7) {x1[]*tc_bright.coef[1] + x2[]*tc_wet.coef[2]+ x3[]*tc_wet.coef[3]+ x4[]*tc_wet.coef[4]+ x5[]*tc_wet.coef[5]+ x7[]*tc_wet.coef[6]},
  filename=paste(out.tcbright,"/","tcbright.",substr(basename(list[i]),1,16),sep=""),overwrite=T)
  overlay(s[[1]],s[[2]],s[[3]],s[[4]],s[[5]],s[[6]],fun=function(x1,x2,x3,x4,x5,x7) {x1[]*tc_green.coef[1] + x2[]*tc_wet.coef[2]+ x3[]*tc_wet.coef[3]+ x4[]*tc_wet.coef[4]+ x5[]*tc_wet.coef[5]+ x7[]*tc_wet.coef[6]},
  filename=paste(out.tcgreen,"/","tcgreen.",substr(basename(list[i]),1,16),sep=""),overwrite=T)
  file.remove(paste(srdir,"/",sr_files,sep=""))
}
