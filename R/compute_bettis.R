#' Calculating Betti numbers of list of point cloud data
#' 
#' @param X list of point cloud data
#' @param maxdim the max dimension to compute persistent homology
#' @param maxscale the maximum of range to compute
#' @param samples the number of subsamples
#' @param const.size the number of points in one subsamples
#'   when 'const.size=F', the number of points in one subsamples is the 80% number of points in original data
#' @export
#' @return list that have 1st and 2nd Betti numbers, the number of data, the number of points in one subsamples
#'   the number of subsamples, max dimension, and max scale

calc_bettis <- function(X, maxdim, maxscale, samples, const.size=F){
  aggr1 <- matrix(0,length(X),1)
  aggr2 <- matrix(0,length(X),1)
  rownames(aggr1) <- paste0("data-set", 1:length(X))
  dimnames(aggr2) <- rownames(aggr1)
  
  for(t in 1:length(X)){
    
    cat("data set", t, "calculating\n")
    if(const.size==F){size<-nrow(X[[t]])*(4/5)}
    else{size<-const.size}
    
    B <- bootstrapper(X[[t]], size, samples)
    speak <- bootstrap_homology(B,maxdim,maxscale)
    m <- sapply(1:maxdim,function(d)speak[[paste0("dim",d,"dhole")]])
    
    aggr1[t,1] <- m[1]
    aggr2[t,1] <- m[2]
  }
  
  aggrs <- list(aggr1,aggr2)
 
  aggrs <- append(aggrs,list(Xsize=sapply(1:length(X), function(l)nrow(X[[l]]),Xsamples=length(X),
                             Bsize=size,Bsamples=samples,
                             maxdim=maxdim,maxscale=maxscale))
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


#'Computing cycles in subsamples
#'

bootstrap_homology <- function(X,maxdim,maxscale,const.band=0,maximum.thresh = F){

  if(!("bootsSamples" %in% class(X))) stop("input must be bootsSamples")
  peak <- matrix(0,maxdim,length(X))
  tseq <- seq(0,maxscale,length.out = 1000)
  diags <- lapply(X,function(x)calcPhom(x,maxdim,maxscale,ret = T,plot = F))
  print(sapply(diags,function(diag)calcDiagCentroid.mk2(diag)[1]))
  band <- ifelse(const.band==0,max(sapply(diags,function(diag)calcDiagCentroid.mk2(diag)[1])),const.band)
  print(band)
  
  for (t in 1:length(X)) {
    land <- lapply(1:maxdim,function(d)landscape(diags[[t]][[1]],dimension = d,KK = 1,tseq = tseq))
    if(maximum.thresh) band <- max(sapply(land,max))/4
    for(d in 1:maxdim){
      peak[d,t] <- calc.landscape.peak(X=land[[d]], thresh = (band/d), tseq=tseq)
    }
  }
  
  dimnames(peak) <- list(paste0("dim",1:maxdim),paste0("sample",1:length(X)))
  bootstrap.summary <- list(peak=peak)
  bootstrap.summary <- append(bootstrap.summary,c(band=band,show.hole.density(peak)))
  class(bootstrap.summary) <- "smoothPhom"
  return(bootstrap.summary)
}