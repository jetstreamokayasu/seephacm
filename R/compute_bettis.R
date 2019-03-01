#' Calculating Betti numbers of list of point cloud data
#' 
#' @param X list of point cloud data
#' @param maxdim the max dimension to compute persistent homology
#' @param maxscale the maximum of range to compute
#' @param samples the number of subsamples
#' @param const.size the number of points in one subsamples
#'   when 'const.size=F', the number of points in one subsamples is the 80% number of points in original data
#' @export
#' @return aggrs list that have 1st and 2nd Betti numbers, the number of data, the number of points in one subsamples
#'   the number of subsamples, max dimension, and max scale

calc_bettis <- function(X, maxdim, maxscale, samples, const.size=F){
  aggr1 <- matrix(0,length(X),1)
  aggr2 <- matrix(0,length(X),1)
  rownames(aggr1) <- paste0("data-set", 1:length(X))
  rownames(aggr2) <- rownames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    if(const.size==F){size<-nrow(X[[t]])*(4/5)}
    else{size<-const.size}
    
    B <- bootstrapper(X[[t]], size, samples)
    speak <- bootstrap_homology(B,maxdim,maxscale)
    
    aggr1[t,1] <- mean(speak[1,])
    aggr2[t,1] <- mean(speak[2,])
  }
  
  aggrs <- list(aggr1,aggr2)
 
  aggrs <- append(aggrs,list(Xsize=sapply(1:length(X), function(l)nrow(X[[l]])), Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim, maxscale=maxscale))
  class(aggrs) <- "bettiComp"
  
  return(aggrs)
}

#' Making subsamples
#' @importFrom TDA hausdInterval
#' @return "bootsSamples" class

bootstrapper <- function(X,size,samples){
  
  hausdint <- TDA::hausdInterval(X = X,m = size,B = samples)
  X <- lapply(1:samples,function(i)X[sample(nrow(X),size),])
  attr(X,"size") <- size
  attr(X,"samples") <- samples
  attr(X,"hausd") <- hausdint
  class(X) <- c("matrix","bootsSamples")
  return(X)
}


#' Computing cycles in subsamples
#' @importFrom phacm compute_pd
#' @importFrom phacm compute_pl
#' @importFrom phacm count_smooth_maximal
#' @importFrom magrittr %>% 

bootstrap_homology <- function(X,maxdim,maxscale,const.band=0,maximum.thresh = F){

  if(!("bootsSamples" %in% class(X))) stop("input must be bootsSamples")
  peak <- matrix(0,maxdim,length(X))

  diags <- lapply(X,function(x)phacm::compute_pd(x,maxdim,maxscale))
  
  for (t in 1:length(X)) {
    pd<-diags[[t]]
    attr(pd, "maxdimension")<-maxdim
    attr(pd, "scale")<-c(0, maxscale)
    land <- phacm::compute_pl(pd)
    if(maximum.thresh) band <- max(sapply(land,max))/4
    peak[, t]<- phacm::count_smooth_maximal(land, exist.method = per_mean, cutoff.method = per_mean)[, "betti"] %>% 
      as.matrix()
  }
  
  dimnames(peak) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))

  return(peak)
}


#' Calculating threshold using the mean of persistence
#' @importFrom phacm zero_hat_double_threshold

per_mean<-function(pd){
  
  mean<-phacm::zero_hat_double_threshold(pd)
  
  return(mean/2)
  
}
