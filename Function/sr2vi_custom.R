sr2vi_custom <- function (x, vi = "ndvi", e = NULL, mask = NULL, keep = c(0), 
                          L = 0.5, ...) 
{
  clean <- function(x, y) {
    x[!(y %in% keep)] <- NA
    return(x)
  }
  if (extension(x[1]) == ".hdf") {
    x <- unlist(sapply(FUN = function(x) {
      try(get_subdatasets(x), silent = TRUE)
    }, X = x), use.names = FALSE)
  }
  if (any(grepl(pattern = sprintf("^.*%s($|\\.tif)", vi), x = x, 
                ignore.case = TRUE))) {
    vi <- raster(grep(pattern = sprintf("^.*%s($|\\.tif)", 
                                        vi), x = x, value = TRUE, ignore.case = TRUE))
    if (!is.null(mask)) {
      mask <- raster(grep(pattern = sprintf("^.*%s($|\\.tif|_band$)", 
                                            mask), x = x, value = TRUE))
    }
    if (!is.null(e)) {
      if (class(e) != "extent") {
        e <- extent(e)
      }
      vi <- crop(vi, e)
      if (!is.null(mask)) {
        mask <- crop(mask, e)
      }
    }
    if (!is.null(mask)) {
      vi <- overlay(x = vi, y = mask, fun = clean, ...)
    }
    return(vi)
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
    ind <- viFormula$ind
    x0 <- grep(pattern = sprintf("^.*(%s)($|\\.tif)", paste(ind, 
                                                            collapse = "|")), x = x, value = TRUE)
    bands <- lapply(X = x0, FUN = raster)
    if (!is.null(mask)) {
      mask <- raster(grep(pattern = sprintf("^.*%s($|\\.tif|_band$)", 
                                            mask), x = x, value = TRUE))
    }
    if (!is.null(e)) {
      if (class(e) != "extent") {
        e <- extent(e)
      }
      bands <- lapply(X = bands, FUN = crop, y = e)
      if (!is.null(mask)) {
        mask <- crop(mask, e)
      }
    }
    fun <- viFormula$fun
    dots <- list(...)
    doListDots <- c(bands, fun = fun, dots)
    doList <- c(bands, fun = fun)
    if (is.null(mask)) {
      vi <- do.call(what = raster::overlay, args = doListDots)
    }
    else {
      previ <- do.call(what = raster::overlay, args = doList)
      vi <- overlay(x = previ, y = mask, fun = clean, ...)
    }
    return(vi)
  }
}