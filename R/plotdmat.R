#' Plot the distance matrix in dbICC 
#' 
#' This function is to plot the distance matrix in dbICC including the segmentation boundary (colored in red) between the within-
#' distance and between-distance based on \code{ggplot2}.
#' 
#' 
#' @param dmat See \code{\link{dm2icc}}.
#' @param nsub number of subject or individual.
#' @param nmea a vector containing number of  measurement for each subject or individual;
#' if \code{nmea} is a scalar, it means each subject shares the same number of the measurement.
#' @param xlab a title for x axis.
#' @param legend logicals. It defaults to show the legend bar.
#' 
#' @return a ggplot object

#' 
#' @author Meng Xu \email{mxu@@campus.haifa.ac.il}
#' 
#' @seealso \code{\link{dm2icc}},\code{\link{dm2icc.bt}}
#' 
#' @references 
#' \itemize{
#' \item Xu, M., Reiss, P. T., and Cribben, I. (2020). Generalized reliability based on distances. Biometrics, to appear. \url{https://arxiv.org/abs/1912.07137}.
#' }
#' @keywords Distance matrix
#' 
#' @import ggplot2 reshape2
#' @export 
#' @examples
#' # See example for dm2icc
#' 
plotdmat <-
function(dmat, nsub, nmea, xlab=NULL, legend=TRUE){
    # plot dist matrix with boundary between within and between distance 
    dmat<-as.matrix(dmat)
    if (length(nmea)==1) nmea=rep(nmea,nsub)
    blkfun<-function(i) sum(nmea[1:i-1])+c(1,nmea[i])
    blkm<-function(i) c(0,-nmea[i]+1)
    xyl<-unlist(lapply(1:nsub,blkfun))
    ym<-unlist(lapply(1:nsub,blkm))
    xyl<-rbind(xyl,xyl+ym)+matrix(rep(c(-.5,-.5,0.5,-0.5),times=nsub),2)
    xyl<-t(xyl)
    colnames(xyl)<-c("x","y")
    
    A<-melt(dmat)
    colnames(A)<-c('x','y','value')
    B<-data.frame(xyl)
    p<-ggplot(A,aes_(~x,~y))+geom_raster(aes_(fill=~value))+
      geom_path(data=B,aes_(~x,~y),color="red")+
      geom_path(data=B,aes_(~y,~x),color="red")+xlab(xlab)+
      theme(panel.grid = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x=element_text(size=18),
            axis.text = element_blank(),
            axis.ticks = element_blank())+scale_y_reverse() 
    if (!legend) p <- p + theme(legend.position="none")
    return(p)
}
