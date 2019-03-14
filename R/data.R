#' Generating points distributing on disk uniformly
#' @param N the number of points
#' @param R radius
#' @importFrom stats runif
#' @export

uniDiskMake<-function(N, R){

  theta <- stats::runif(N, min=0, max=2*pi)
  r <- sqrt(2*runif(N, min=0, max=0.5*R^2))
  df <- data.frame(x=r*cos(theta), y=r*sin(theta))
  return(df)

}


#' Combing 'plot3d' and 'aspect3d' to plot 3d data
#' @param X 3d data
#' @importFrom rgl plot3d
#' @importFrom rgl aspect3d
#' @export

figurePlot3d<-function(X){

  rgl::plot3d(X)
  rgl::aspect3d("iso")

}

#' Saving variable into RData file
#' @param ... one variable
#' @importFrom magrittr %>%
#' @export

save2Rdata <- function(...) {

  elp <- list(...)
  elname <- substitute(...) %>% as.character()
  assign(elname, elp[[1]])
  save(list = elname, file = paste0("./data/", gsub("\\.", "_", elname), ".RData"))
}
