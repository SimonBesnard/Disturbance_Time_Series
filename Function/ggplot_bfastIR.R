ggplot.bfastIR <- function(x, seg = TRUE, order, formula, numbering = FALSE) {
  
  ggdf <- x$df
  ggdf[,'breaks'] <- NA
  ggdf$breaks[x$breaks$breakpoints] <- 1
  
  xIntercept <- ggdf$time[ggdf$breaks == 1]
  ggdf[,'breakNumber'] <- NA
  if (!is.na(x$breaks$breakpoints)) {
    ggdf$breakNumber[!is.na(ggdf$breaks)] <- 1:length(x$breaks$breakpoints)
  } else {
    ggdf$breakNumber[floor(nrow(ggdf)/2)] <- "No Break"
  }
  ggdf[,'maxY'] <- max(ggdf$response)
  
  
  gg <- ggplot(ggdf, aes(time, response)) +
    geom_line(size=0.2) +
    geom_point(color = 'green', size=1 ) +
    geom_vline(xintercept = xIntercept, color = 'red', linetype = 'dashed')
  if(numbering) {
    gg <- gg + geom_text(aes(x = time + 0.5, y = maxY, label = breakNumber))
  }
  
  gg <- gg +
    scale_x_continuous(breaks=floor(min(ggdf$time)):ceiling(max(ggdf$time))) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  
  if(seg && !is.na(x$breaks$breakpoints)) {
    # Segments on time column
    segments <- c(ggdf$time[c(1,x$breaks$breakpoints, nrow(ggdf))])
    for(i in seq_along(segments[-1])) {
      predTs <- bfastts(rep(NA, ncol(ggdf)), date_decimal(ggdf$time), type = 'irregular')
      predDf <- bfastpp(predTs, order = order, na.action = na.pass)
      predDfSub <- subset(predDf, time <= segments[i + 1] & time >= segments[i])
      trainDfSub <- subset(ggdf, time <= segments[i + 1] & time >= segments[i])
      model <- lm(formula = formula, data = trainDfSub)
      predDfSub$pred <- predict(model, newdata = predDfSub)
      
      gg <- gg + geom_line(data = predDfSub, aes(x = time, y = pred), color = 'blue', size=0.2)
      
      
    }    
  }
  gg
}