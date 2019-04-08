#' Distance-based Intraclass Correlation Coefficient (dbICC)
#' 
#' When the distance among any test-retest datatype can be defined, its reliability 
#' performed in Reiss et al. (2019) can be obtained by this function, to calculate the 
#' value of dbICC through the related distance matrix. 
#' 
#' 
#' 
#' @param dmat A distance matrix or a object of \code{dist}, of dimension \code{sum(nmea)*sum(nmea)}. Note that the order of the rows or
#' columns in the distance matrix is grouped of each subject or individual. The details refer to Figure 1 of
#' Reiss el at., 2019.
#' @param nsub Number of subject or individual.
#' @param nmea A vector containing number of measurement for each subject or individual; 
#' if \code{nmea} is a scalar, it means each subject shares the same number of the measurement.
#' 
#' @return A scalar, giving the dbICC value
#' 
#' @author Meng Xu \email{mxu@@campus.haifa.ac.il}, Philip Reiss
#' 
#' @seealso \code{\link{plotdmat}},\code{\link{boot.dbicc}}
#' 
#' @references
#' \itemize{
#' \item P.T. Reiss, M. Xu, I. Cribben (2019). Generalized reliability based on distances. To be submitted.
#' }
#' 
#' @keywords dbICC, reliability
#' @import Matrix MASS
#' @export
#' @examples
#' 
#' ##Point estimates of dbICC
#' 
#' # Generation function for R^2 points from multi-normal distribution
#' R2gen<-function(nsub,nmea,m=1,variance,sds=NULL,pt=FALSE){
#'     if (is.null(sds)==FALSE) set.seed(sds)
#'     if (length(nmea)==1) nmea<-rep(1,nsub)*nmea
#'     sig1<-diag(rep(variance,2))
#'     mu<-c(0,0)
#'     sig2<-diag(rep(1,2)) # true-value variance of the 2-d normal distribution
#'     t<-MASS::mvrnorm(nsub,mu,sig2)
#'     e<-MASS::mvrnorm(sum(nmea),mu,sig1/m)
#'     p<-matrix(apply(t,2,rep,times=nmea),ncol=2)+e#(I*J)x2
#'     if (pt==TRUE) return(list(t=t,p=p))
#'     if (pt==FALSE) return(p)
#' }
#' 
#' # set the number of the point
#' I <- 10
#' 
#' # set the number of the measurement for each point
#' J <- 4
#' 
#' # generate the sample of R^2 points
#' varl <- .25 # error variance of the 2-d normal distribution
#' pij <- R2gen(I,J,variance=varl)
#' 
#' # calculate the squared distance matrix via Euclidean distance
#' distmat<-as.matrix(dist(pij))
#' 
#' #plot the distance matrix
#' plotdmat(distmat,I,J)
#' 
#' # dbICC value
#' dm2icc(distmat,I,J)
#' 
#' 
 
dm2icc <-
function(dmat, nsub, nmea) {
  #nmea: the number of the measurement
  if (length(nmea)==1) nmea<-rep(nmea,nsub)
  dmat<-as.matrix(dmat^2)
  wmask <- as.matrix(Matrix::bdiag(lapply(nmea, function(i) matrix(1, i, i))))
  diag(wmask) <- NA
  bmask <- matrix(1, sum(nmea), sum(nmea)) - wmask
  wmask[wmask == 0] <- bmask[bmask == 0] <- NA
  bmask.ah <- bmask
  dboot<-dmat
  1 - mean(dboot * wmask, na.rm = TRUE) / mean(dboot * bmask.ah, na.rm = TRUE)
}
