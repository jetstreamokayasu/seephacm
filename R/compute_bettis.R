#' Calculating Betti numbers of list of point cloud data
#'
#' @param X list of point cloud data
#' @param maxdim the max dimension to compute persistent homology
#' @param maxscale the maximum of range to compute
#' @param samples the number of subsamples
#' @param const.size the number of points in one subsamples.
#'   when 'const.size=F', the number of points in one subsamples is the 80% number of points in original data
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
#' @param X origin data
#' @param size the number of data points in one subsample
#' @param samples the number of subsamples
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
#' @param n dimension
#' @importFrom phacm persistence
#' @export
calcper<-function(pd, n){

  per<-phacm::persistence(pd)
  per.dim<-per[per[,"dim"]==n, "persistence"] %>% as.vector()

  if(n==0){per.dim<-per.dim[-1, ]}

  return(per.dim)
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
#' @param line draw threshold in a graph
#' @importFrom TDA landscape
#' @importFrom graphics plot
#' @export

calc_landscape<-function(diag, maxscale, line=T){

  thresh<-phacm::calcDiagCentroid(diag)[3]
  tseq <- seq(0, maxscale, length = 1000)
  Land.dim1 <- TDA::landscape(diag[[1]], dimension = 1, KK = 1, tseq)
  graphics::plot(tseq, Land.dim1, type = "l", col=2, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(Land.dim1)+1)/2), main ="1-degree landscape")
  if(line){
    graphics::abline(h=thresh)
  }

  if(length(diag[[1]][diag[[1]][,1]==2,])>0){
    Land.dim2 <- TDA::landscape(diag[[1]], dimension = 2, KK = 1, tseq)
    graphics::plot(tseq, Land.dim2, type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(Land.dim2)+1)/2), main ="2-degree landscape")
    if(line){
      graphics::abline(h=thresh/2)
    }

    return(list(tseq=tseq, Land.dim1=Land.dim1, Land.dim2=Land.dim2, thresh=thresh))

  }else{return(list(tseq=tseq, Land.dim1=Land.dim1, thresh=thresh))}
}


#' Plotingpersistence  landscape
#' @param land persistence landscape calculated by 'calc_landscape'
#' @export

plotLandscape<-function(land){

  plotland<-lapply(2:(length(land)-1), function(k){

    if(names(land)[k]=="Land.dim1"){

      graphics::plot(land[[1]], land[[k]], type = "l", col=k, xlab = "(Birth + Death) / 2", ylab = "(Death - Birth) / 2", ylim=c(0, round(max(land[[k]])+1)/2), main =paste0(k-1, "-degree landscape"))
      graphics::abline(h=land[["thresh"]])

    }

    else if(names(land)[k]=="Land.dim2"){

      graphics::plot(land[[1]], land[[k]], type = "l", col=3, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(land[[k]])+1)/2), main =paste0(2, "-degree landscape"))
      graphics::abline(h=land[["thresh"]]/2)

    }else{

      graphics::plot(land[[1]], land[[k]], type = "l", col=k, xlab = "(Birth + Death) / 2",ylab = "(Death - Birth) / 2", ylim=c(0, round(max(land[[k]])+1)/2), main =paste0(k-1, "-degree landscape"))

    }

  })

}

