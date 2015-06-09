# NDVI --------------------------------------------------------------------- 
.ndvi <- function() { 
  ind <- c('band3','band4') 
  fun <- function(x1, x2) { 
    ndvi <- 10000 * (x2 - x1)/(x2 + x1) 
    return(ndvi) 
  } 
  return(list(ind=ind,fun=fun)) 
} 
18 

19 

# EVI --------------------------------------------------------------------- 
.evi <- function() {
  ind <- c('band1','band3','band4') 
  fun <- function(x1, x3, x4){  
    evi <- 10000 * 2.5 * (x4/10000 - x3/10000)/(x4/10000 + 6 * x3/10000 - 7.5 * x1/10000 + 1) 
    return(evi) 
  } 
  return(list(ind=ind, fun=fun)) 
} 

# NBR --------------------------------------------------------------------- 
.nbr <- function() { 
  ind <- c('band4','band7') 
  fun <- function(x1, x2) { 
    ndvi <- 10000 * (x1 - x2)/(x1 + x2) 
    return(ndvi) 
  } 
  return(list(ind=ind, fun=fun)) 
} 


# Tasseled Cap Components -------------------------------------------------  

.tcbright <- function() {  
  ind <- c('band1','band2','band3','band4','band5','band7')   
  # make compatible with getSceneinfo() output   
  tc_coef <- c(0.2043, 0.4158, 0.5524, 0.5741, 0.3124, 0.2303)  
  fun <- function(x1, x2, x3, x4, x5, x7) {  
    tcbright <- sum(c(x1, x2, x3, x4, x5, x7) * tc_coef)  
  }  
  return(list(ind=ind,fun=fun)) 
}  

.tcgreen <- function() {  
  ind <- c('band1','band2','band3','band4','band5','band7')  
  # make compatible with getSceneinfo() output  
  tc_coef <- c(-0.1603, -0.2819, -0.4934, 0.7940, -0.0002, -0.1446)  
  
  fun <- function(x1, x2, x3, x4, x5, x7) {  
    tcgreen <- sum(c(x1, x2, x3, x4, x5, x7) * tc_coef)  
  }  
  return(list(ind=ind,fun=fun)) 
}  

.tcwet <- function() {  
  ind <- c('band1','band2','band3','band4','band5','band7')  
  # make compatible with getSceneinfo() output  
  tc_coef <- c(0.0315,  0.2021,  0.3102,  0.1594, -0.6806, -0.6109)  
  
  fun <- function(x1, x2, x3, x4, x5, x7) {  
    tcwet <- sum(c(x1, x2, x3, x4, x5, x7) * tc_coef)  
  }  
  return(list(ind=ind,fun=fun))
}  
