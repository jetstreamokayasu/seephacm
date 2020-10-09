#' Generating points distributing on disk uniformly
#' @param N the number of points
#' @param R radius
#' @importFrom stats runif
#' @export

diskUnif<-function(N, R){

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

#' Saving variable into RData file in "data" directory under working directory.
#' When there isn't "data" directory in working directory, this function make "data" directory.
#' @param ... one variable
#' @importFrom magrittr %>%
#' @export

save2RData <- function(...) {

  elp <- list(...)
  elname <- substitute(...) %>% as.character()
  assign(elname, elp[[1]])

  dir<-match("data", list.files())
  if(is.na(dir)){
    dir.create("./data")
  }else{
    if(!(file.info("./data")$isdir)){dir.create("./data")}
  }

  save(list = elname, file = paste0("./data/", elname, ".RData"))

}


#'This function samples n points from the ellipsoid with equation (x/a)^2+(y/b)^2+(z/c)^2=1, uniformly with respect to its surface.
#'This function is tentative. There is no guarantee that points are sampled unifomly.
#'@param n the number of points in the sample.
#'@param a the diameter of x-axis direction.
#'@param b the diameter of y-axis direction.
#'@param c the diameter of z-axis direction.
#'@importFrom stats runif
#'@export
#'

xEllip_unif<-function(n, a, b, c){

  p<-stats::runif(n, -1, 1)
  theta<-acos(-p)
  phi<-stats::runif(n, 0, 2*pi)

  x<-a*sin(theta)*cos(phi)
  y<-b*sin(theta)*sin(phi)
  z<-c*cos(theta)

  return(cbind(x, y, z))

}

#'The function for calculate degrees of the uniform distribution respect to 3D-torus.
#'This function is tentative.
#'@param r the smallest radius of 3D torus.
#'@param R1 the largest radius of 3D torus.
#'@param R2 the medim radius of 3D torus.
#'@param theta an angle of 3D torus.
#'@param phi an angle of 3D torus.
#'
rad_distribute<-function(r, R1, R2, theta, phi){

  g<-1/R1*(1 + (r/R2)*sin(theta))*(1/R2 + (1/R1)*sin(phi)*(1 + (r/R2)*sin(theta)))

  return(g)

}

#'This function samples n points from the surface or the 3D torus uniformly.
#'This function is tentative. There is no guarantee that points are sampled unifomly.
#'@param n the number of points in the sample.
#'@param r the smallest radius of 3D torus.
#'@param R1 the largest radius of 3D torus.
#'@param R2 the medim radius of 3D torus.
#'@importFrom stats runif
#'@export

x3Dtorus_unif<-function(n, r, R1, R2){

  theta<-c()
  phi<-c()

  tp<-expand.grid(seq(0, 2*pi, length=100), seq(0, 2*pi, length=100))

  g<-apply(tp, 1, function(i){

    g0<-rad_distribute(r, R1, R2, i[1], i[2])

    return(g0)

  })

  max_g<-rad_distribute(r, R1, R2, pi/2, pi/2)
  min_g<-min(g)

  while(length(theta) < n){

    xvec<-stats::runif(1, 0, 2*pi)
    yvec<-stats::runif(1, 0, 2*pi)
    zvec<-stats::runif(1, min_g, max_g)

    fxy<-rad_distribute(r, R1, R2, xvec, yvec)

    if(zvec <= fxy){
      theta<-c(theta, xvec)
      phi<-c(phi, yvec)
    }

  }

  delta<-stats::runif(n, 0, 2*pi)

  x<-sin(delta)*(R1 + sin(phi)*(R2 + r*sin(theta)))
  y<-cos(delta)*(R1 + sin(phi)*(R2 + r*sin(theta)))
  z<-cos(phi)*(R2 + r*sin(theta))
  w<-r*cos(theta)

  out<-cbind(x, y, z, w)

  return(out)

}
