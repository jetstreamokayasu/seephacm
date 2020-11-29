#'Plotting persistence diagram in a list
#'
#'@param diags list of persisitence diagram
#'@importFrom graphics plot
#'@importFrom graphics par
#'@export
#'
plot_diagrams<-function(diags){

  graphics::par(mfrow=c(((length(diags)%/%4)+1), 4))
  for (k in 1:length(diags)) {
    graphics::plot(diags[[k]][[1]], cex.lab=1.6, cex.axis=1.6)
  }

  graphics::par(mfrow=c(1, 1))
}

#'Plotting persistence landscape in a list
#'
#'@param lands list of persistence landscape
#'@param dim dimension to plot
#'@param xlim the x limits (x1, x2) of the plot of persistence landscape
#'@param ylim the y limits of the plot of persistence landscape
#'@importFrom graphics par
#'@export
#'
plot_lands<-function(lands, dim, xlim, ylim){

oldpar<-graphics::par(no.readonly=T)
graphics::par(cex.lab=2, cex.main=2, cex.axis=2, plt = c(0.2, 0.9, 0.2, 0.9))

graphics::par(mfrow=c(2, (length(lands)%/%2+length(lands)%%2)))

for (k in 1:length(lands)) {
  plot_landscape(land = lands[[k]], dim = dim, xlim = xlim, ylim = ylim)
}

graphics::par(oldpar)

}

#'Plotting persistence landscape
#'
#'@param land persistence landscape calculated by 'calc_landscape'
#'@param dim dimension of persistence landscape to plot
#'@param xlim the x limits (x1, x2) of the plot of persistence landscape
#'@param ylim the y limits of the plot of persistence landscape
#'@importFrom graphics plot
#'@importFrom graphics abline
#'@export
#'

plot_landscape<-function(land, dim, xlim, ylim){

  if(missing(xlim)){xlim<-c(0, max(land[["tseq"]]) )}
  if(missing(ylim)){ylim<-c(0, round(max(land[[paste0(dim, "-land")]])+0.05, digits = 1))}

  graphics::plot(land[["tseq"]], land[[paste0(dim, "-land")]], type = "l", col=dim+1, xlab = "(Birth + Death) / 2", ylab = "(Death - Birth) / 2", xlim=xlim, ylim=ylim, main = paste0(dim, "-degree landscape"))
  graphics::abline(h=land[["thresh"]]*((2*pi)/surface_nshpere(dim)))


}

#'Draw cicle in a figure with 'polygon'.
#'@param x x-coordinate of the center of a circle.
#'@param y y-coordinate of the center of a circle.
#'@param r the radius of a circle.
#'@param col color.
#'@importFrom graphics polygon
#'@importFrom grDevices rgb
#'@export
#'

plot_circle<-function(x, y, r, col){

  if(missing(col)){col=grDevices::rgb(1, 1, 1, 0)}

  theta <- seq(-pi, pi, length=100)
  graphics::polygon(x + r*cos(theta), y + r*sin(theta), col=col)


}


