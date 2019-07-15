#'Plotting persistence diagram in a list
#'
#'@param diags list of persisitence diagram
#'@importFrom graphics plot
#'@importFrom graphics par
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
#'@importFrom graphics par
#'
plot_lands<-function(lands, dim){

oldpar<-graphics::par(no.readonly=T)
graphics::par(cex.lab=2, cex.main=2, cex.axis=2, plt = c(0.2, 0.9, 0.2, 0.9))

graphics::par(mfrow=c(2, (length(lands)%/%2+length(lands)%%2)))

for (k in 1:length(lands)) {
  plotLandscape(lands[[k]][dim-4])
}

graphics::par(oldpar)

}
