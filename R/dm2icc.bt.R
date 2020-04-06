#' Bootstrap Confidence intervals for dbICC
#' 
#' Nonparametric bootstrapping can be used to construct confidence 
#' intervals for the Distance-based Intraclass Correlation Coefficient 
#' (dbICC) based on samples of the subjects with replacement.
#' 
#' @param dmat a distance matrix or an object of \code{dist}, of dimension \code{sum(nmea)*sum(nmea)}. 
#' Note that the structure of the distance matrix, with the rows or
#' columns is grouped by subjects or individuals. The details refer to Figure 4 and Table 2 of
#' Xu el at.(2020).
#' @param nsub number of subject or individual.
#' @param nmea a vector containing number of the measurement for each subject or individual; 
#' if \code{nmea} is a scalar, it means each subject shares the same number of the measurement.
#' @param nB number of bootstrap.
#' 
#' @param adhoc a logical variable, whether to apply the ad hoc correction when
#' estimating the dbICC from a bootstrap sample. Default is \code{TRUE}.
#' @param probs a vector of probabilities with values in [0,1]. \code{c(.025, .975)} is default.
#' @return estimates of underlying dbicc sample quantiles in \code{probs}. 95% (\code{c(.025, .975)}) is default.
#' 
#' @author Meng Xu \email{mxu@@campus.haifa.ac.il}
#' 
#' @seealso \code{\link{plotdmat}},\code{\link{dm2icc}}
#' 
#' @references
#' \itemize{
#' \item Xu, M., Reiss, P. T., and Cribben, I. (2020). Generalized reliability based on distances. Biometrics, to appear. \url{https://arxiv.org/abs/1912.07137}.
#' }
#' 
#' @keywords Bootstrap, dbICC, reliability
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
#' distmat<-as.matrix(dist(pij))
#' 
#' ##Bootstrap Confident intervals
#' 
#' dm2icc.bt(distmat, I, J, 500)


dm2icc.bt <-
  function(dmat, nsub, nmea, nB, probs=c(.025, .975), adhoc = TRUE){
    bicc <- numeric(nB)
    for (b in 1:nB){
      bs <- sort(sample(nsub, replace = TRUE))
      bicc[b] <- boot.dbicc(dmat, nsub, nmea, bootsamp=bs, adhoc=adhoc)
    }
    return(stats::quantile(bicc, probs=probs))
  }