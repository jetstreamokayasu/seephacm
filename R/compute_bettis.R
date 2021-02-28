#' Calculating Betti numbers of list of point cloud data
#'
#' @param X list of point cloud data
#' @param maxdim the max dimension to compute persistent homology
#' @param maxscale the maximum of range to compute
#' @param samples the number of subsamples
#' @param const.size the number of points in one subsamples.
#'   when 'const.size=F', the number of points in one subsamples is the 80\% number of points in original data
#' @export
#' @return return list that have 1st and 2nd Betti numbers, the number of data, the number of points in one subsamples
#'   the number of subsamples, max dimension, and max scale

calc_bettis <- function(X, maxdim, maxscale, samples, const.size=F){
  aggr1 <- matrix(0,length(X),1)
  aggr2 <- matrix(0,length(X),1)
  rownames(aggr1) <- paste0("data-set", 1:length(X))
  rownames(aggr2) <- rownames(aggr1)

  for(t in 1:length(X)){

    cat("data set", t, "calculating\n")
    if(const.size==F){size<-nrow(X[[t]][[2]])*(4/5)}
    else{size<-const.size}

    B <- bootstrapper(X[[t]][[2]], size, samples)
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
#' @param X origin data
#' @param size the number of data points in one subsample
#' @param samples the number of subsamples
#' @return "bootsSamples" class

bootstrapper <- function(X,size,samples){

  hausdint <- TDA::hausdInterval(X = X, m = size, B = samples)
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
#' @param X subsamples
#' @param maxdim the max dimension to compute persistent homology
#' @param maxscale the maximum of range to compute
#' @param const.band constant band to discriminate cycles and noises
#' @param maximum.thresh use a quarter of max persistence as threshold

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
#' @param pd persisitence diagram by 'phacm'
#' @export

per_mean<-function(pd){

  mean<-phacm::zero_hat_double_threshold(pd)

  return(mean/2)

}


#' Calculating persistence of n-dimensional holes
#' @param pd persisitence diagram by 'phacm'
#' @param dim dimension
#' @importFrom assertthat assert_that
#' @export
calc_per<-function (pd, dim){

  assertthat::assert_that((length(dim) == 1) && is.numeric(dim))

  pers <- pd[pd[, 1] == dim, 3] - pd[pd[, 1] == dim, 2]

  attr(pers, "pers_dim") <- dim

  return(pers)

}


#' Counting the number of data sets estimated that have n holes
#' @param aggr estimated result using 'calc_bettis'
#' @param dim dimension of holes
#' @export

cycle_number<-function(aggr, dim){

  cyc.num<-max(aggr[[dim]])

  if(cyc.num > 4){col<-cyc.num}
  else{col<-4}

  estimate<-matrix(0, ncol(aggr[[dim]]), col+1)
  colnames(estimate)<-c(paste0(0:col, "hole"))

  for (i in 1:(col+1)) {
    for (j in 1:ncol(aggr[[dim]])) {
      estimate[j, i]<-meetNumber(aggr[[dim]][, j], i-1)
    }
  }

  return(estimate)

}


#' Counting correct result
#' @param result estimated Betti numbers
#' @param correct correct Betti numbers

meetNumber<-function(result, correct){

  return(length(result[result >= (correct-0.5) & result < (correct+0.5)]))

}

#' Calculating succsess rate of results
#' @param aggrlist list of estimate result of 'calc_bettis'
#' @param correct vector of correct Betti numbers
#' @export
#'

aggr_success_rates<-function(aggrlist, correct){

  rates<-lapply(aggrlist, function(aggr){

    set.sum<-length(aggr[[1]])

    cycle1<-cycle_number(aggr, 1)
    rate1<-cycle1[correct[1]+1]/set.sum

    cycle2<-cycle_number(aggr, 2)
    rate2<-cycle2[correct[2]+1]/set.sum

    return(list(dim1rate=rate1, dim2rate=rate2))

  })


  return(rates)

}

#' Calculating persistence landscape from persistence diagram by 'TDA' and ploting
#' @param diag persistence diagram
#' @param maxscale the maximum of range to compute
#' @param plot default TRUE; if TRUE, plot persistence landscape
#' @param line default TRUE; if TRUE, draw threshold in a graph
#' @importFrom TDA landscape
#' @importFrom graphics plot
#' @export
#' @return persistence diagram

calc_landscape<-function(diag, maxscale, plot=T, line=T){

  if(is.list(diag)){diag<-diag[["diagram"]]}

  thresh<-persistence_weighted_mean(diag)
  tseq <- seq(0, maxscale, length = 1000) #domain
  maxdim<-max(diag[,1])
  if(missing(maxscale)){maxscale<-max(diag[,3])}

  lands<-lapply(1:maxdim, function(k){

    land<-landscape(diag, dimension = k, KK = 1, tseq)

    if(plot){

      plot(tseq, land, type = "l", col=k+1, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(land)+0.5)), main = paste0(k, "-degree landscape"))
      if(line){abline(h=thresh*((2*pi)/surface_nshpere(k)))}

    }

    attr(land, "maxscale")<-maxscale
    return(land)

  })

  names(lands)<-sapply(1:maxdim, function(k)paste0(k, "-land"))

  return(append(lands, list(thresh=thresh, tseq=tseq)))

}

#'Calculating the mean of persistence except 0-persistence.
#'@param diag a persistence diagram.
#'@importFrom magrittr %>%
#'@return the mean of persistence.
#'

persistence_mean<-function(diag){

  if(class(diag)=="list"){diag<-diag[["diagram"]]}

  maxdim<-max(diag[,1])
  diag <- diag[-which(diag[,1]==0),]

  centroid<-diag[,3]-diag[,2] %>% mean()

  return(centroid)

}

#'Calculate the weighted mean of persistence.
#'n-dimensional persistence multiplied by (the surface area of n-unit sphere/2*pi).
#'@param diag persistence diagram.
#'@importFrom magrittr %>%
#'@return the weighted mean of persistence.
#'

persistence_weighted_mean<-function(diag){

  if(is.list(diag)){diag<-diag[["diagram"]]}

  maxdim<-max(diag[,1])
  diag <- diag[-which(diag[,1]==0),]

  centroid<-lapply(1:maxdim, function(k){

    per<-(diag[diag[,1]==k,3]-diag[diag[,1]==k,2])*(surface_nshpere(k)/(2*pi))
    return(per)

  })

  cpers<-unlist(centroid) %>% mean()

  names(cpers)<-"cpersistence"

  return(cpers)

}


#'Calculate surface area of a n-dimensional sphere
#'@param n dimension of the sphere
#'@param r radius of the sphere
#'@return surface area of a n-dimensional sphere
#'

surface_nshpere<-function(n, r=1){

  s<-(2*pi^((n+1)/2)*(r^n))/gamma((n+1)/2)

  return(s)

}
