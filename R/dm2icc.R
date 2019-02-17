#' Distance-based Intraclass Correlation Coefficient (dbICC)
#' 
#' When the distance among any test-retest datatype can be defined, its reliability 
#' performed in Reiss et al. (2019) can be obtained by this function, to calculate the 
#' value of dbICC through the related squared distance matrix. A sample of the certain subjects with 
#' the replacement, which is used in the boostrap confidence intervals, can be reached as well.  
#' 
#' 
#' @param dmat A squared distance matrix, of dimension \code{sum(nmea)*sum(nmea)}. Note that the order of the rows or
#' columns in the distance matrix is sorted by groups of measurements for each subject or individual. The details refer to
#' Reiss el at., 2019.
#' @param nsub Number of the subject or individual.
#' @param nmea A vector containing number of the measurement for each subject or individual; 
#' if \code{nmea} is a scalar, it means each subject shares the same number of the measurement.
#' @param bootsamp option to a vector requiring a sample with the replacement. Default 
#' includes all the elements of the specified subjects or individuals.
#' @param adhoc option to a logical variable: whether to apply the ad hoc correction when
#' estimating the dbICC from a bootstrap sample. Default is \code{FALSE}.
#' @return A scalar, giving the dbICC value
#' @note %% ~~further notes~~
#' 
#' @author Philip Reiss \email{reiss@@stat.haifa.ac.il}, Meng Xu \email{mxu@@campus.haifa.ac.il}
#' 
#' @seealso \code{\link{plotdmat}}
#' 
#' @references
#' \itemize{
#' \item Reiss et al. (2019). Generalized test-retest reliability based on distances. 
#' }
#' 
#' @keywords dbICC, distance matrix, reliability
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
#'     sig2<-diag(rep(1,2))
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
#' varl <- .25 # variance of the 2-d normal distribution
#' pij <- R2gen(I,J,variance=varl)
#' 
#' # calculate the squared distance matrix via Euclidean distance
#' distmat<-as.matrix(dist(pij)^2)
#' 
#' #plot the squared distance matrix
#' plotdmat(distmat,I,J)
#' 
#' # dbICC value
#' dm2icc(distmat,I,J)
#' 
#' ##Bootstrap Confident intervals
#' 
#' B <- 500
#' bicc <- numeric(B)
#' for (b in 1:B){
#'     bs <- sort(sample(I, replace = TRUE))
#'     bicc[b] <- dm2icc(distmat,I,J,bootsamp=bs)
#' }
#' 
#' quantile(bicc, c(.025, .975))
#' 
 
dm2icc <-
function(dmat, nsub, nmea, bootsamp = 1:nsub, adhoc = FALSE) {
  #nmea: the number of the measurement
  if (length(nmea)==1) nmea<-rep(nmea,nsub)
  wmask <- as.matrix(Matrix::bdiag(lapply(nmea, function(i) matrix(1, i, i))))
  diag(wmask) <- NA
  bmask <- matrix(1, sum(nmea), sum(nmea)) - wmask
  wmask[wmask == 0] <- bmask[bmask == 0] <- NA
  bmask.ah <- bmask
  if (!identical(bootsamp, 1:nsub)){
    wibs <- sapply(bootsamp, function(i) sum(nmea[1:i - 1]) + 1:nmea[i])
    dboot <- dmat[wibs, wibs]
    if (adhoc) {
      for (ii in 1:(nsub - 1))
        for (jj in (ii + 1):nsub) {
          if (bootsamp[ii] == bootsamp[jj]) bmask.ah[(sum(nmea[1:ii - 1]) + 1:nmea[ii]), (sum(nmea[1:jj - 1]) + 1:nmea[jj])] <- bmask.ah[(sum(nmea[1:jj - 1]) + 1:nmea[jj]), (sum(nmea[1:ii - 1]) + 1:nmea[ii])] <- NA
        }
    }
  }
  else dboot<-dmat
  1 - mean(dboot * wmask, na.rm = TRUE) / mean(dboot * bmask.ah, na.rm = TRUE)
}
