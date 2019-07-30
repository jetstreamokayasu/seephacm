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
#'@importFrom graphics par
#'@export
#'
plot_lands<-function(lands, dim){

oldpar<-graphics::par(no.readonly=T)
graphics::par(cex.lab=2, cex.main=2, cex.axis=2, plt = c(0.2, 0.9, 0.2, 0.9))

graphics::par(mfrow=c(2, (length(lands)%/%2+length(lands)%%2)))

for (k in 1:length(lands)) {
  plot_landscape(lands[[k]][dim-4])
}

graphics::par(oldpar)

}

#'Plotting persistence landscape
#'
#'@param land persistence landscape
#'@importFrom graphics plot
#'@importFrom graphics abline
#'@export
#'
plot_landscape<-function(land){

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
