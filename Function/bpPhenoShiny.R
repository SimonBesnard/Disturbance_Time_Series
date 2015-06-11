bpPhenoShiny <- function(x, order=2, formula = response ~ trend + harmon, breaks = NULL, h = 0.15) {
  
  # def function bp2df
  bp2df <- function(bpOut, pp) {
    if(class(bpOut) == 'try-error') {
      return(NULL)
    }
    if(is.na(bpOut$breakpoints)) {
      return(NULL)
    }
    
    nBreak <- length(bpOut$breakpoints)
    
    mat <- matrix(nrow = nBreak, ncol = 11)
    df <- as.data.frame(mat)
    colnames(df) <- c('Date', 'SegDuration_pre', 'SegSlope_pre', 'SegYBegin_pre', 'SegYEnd_pre', 'SegMean_pre', 'SegDuration_post', 'SegSlope_post', 'SegYBegin_post', 'SegYEnd_post', 'SegMean_post')
    
    # Input breaks timing into output df
    df[,1] <- pp$time[bpOut$breakpoints]
    # Input Segments properties
    segments <- c(pp$time[c(1,bpOut$breakpoints, nrow(pp))])
    for (i in 1:(length(segments) - 1)) {
      # To deal with begining of time-series boundary
      if (i == 1) {
        subDf <- subset(pp, time <= segments[i + 1])
      } else {
        subDf <- subset(pp, time <= segments[i + 1] & time > segments[i]) 
      }
      m <- mean(subDf$response, na.rm = TRUE)
      model <- lm(response ~ trend, data = subDf)
      sl <- model$coefficients[2]
      n <- nrow(subDf)
      yBegin <- predict(model, newdata = subDf[1,]) # Beginning of segment (!= to intersection with break)
      yEnd <- predict(model, newdata = subDf[n,]) # End of segment (== intersection with breakLine)
      duration <- diff(range(subDf$time))
      # Fill dataframe 
      if (i > 1) {
        df[i-1, c(7:11)] <- c(sl, yBegin, yEnd, duration, m)
      }
      if (i < (length(segments) - 1))
        df[i, c(2:6)] <- c(sl, yBegin, yEnd, duration, m)
    }
    return(df)        
  }
  
  # def function bp
  bp <- function(ts, order, formula, breaks, h) {
    pp <- bfastpp(ts, order = order)
    bpOut <- try(breakpoints(formula = formula, data = pp, breaks = breaks, h = h))
    dfOut <- bp2df(bpOut, pp) # Returned variable
    out <- list(df = pp,
                breaks = bpOut,
                dfOut = dfOut)
    class(out) <- 'bfastIR'
    return(out)
  }
  
  # def function wrapping the two above functions
  
  ts <- bfastts(data = x, dates = index(x), type = 'irregular')
  out <- bp(ts = ts, order = order, formula = formula, breaks = breaks, h = h)
  return(out)
}
