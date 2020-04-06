#' Bootstrap Confidence intervals for dbICC (internal)
#' 
#' Nonparametric bootstrapping can be used to construct confidence 
#' intervals for the Distance-based Intraclass Correlation Coefficient 
#' (dbICC) based on samples of the subjects with replacement.
#' 
#' @param dmat A distance matrix or an object of \code{dist}, of dimension \code{sum(nmea)*sum(nmea)}. 
#' Note that the structure of the distance matrix, with the rows or
#' columns is grouped by subjects or individuals. The details refer to Figure 1 of
#' Xu el at., 2020.
#' @param nsub Number of subject or individual.
#' @param nmea A vector containing number of the measurement for each subject or individual; 
#' if \code{nmea} is a scalar, it means each subject shares the same number of the measurement.
#' @param bootsamp A sample with replacement.
#' 
#' @param adhoc A logical variable, whether to apply the ad hoc correction when
#' estimating the dbICC from a bootstrap sample. Default is \code{FALSE}.
#' @return A scalar, giving the dbICC value
#' 
#' @author Meng Xu \email{mxu@@campus.haifa.ac.il}, Philip Reiss
#' 
#' @seealso \code{\link{plotdmat}},\code{\link{dm2icc}}
#' 
#' @keywords internal
#' @import Matrix MASS

 
boot.dbicc <-
function(dmat, nsub, nmea, bootsamp, adhoc = FALSE) {
  #nmea: the number of the measurement
  if (length(nmea)==1) nmea<-rep(nmea,nsub)
  dmat<-as.matrix(dmat^2)
  wmask <- as.matrix(Matrix::bdiag(lapply(nmea, function(i) matrix(1, i, i))))
  diag(wmask) <- NA
  bmask <- matrix(1, sum(nmea), sum(nmea)) - wmask
  wmask[wmask == 0] <- bmask[bmask == 0] <- NA
  bmask.ah <- bmask
  if (!all(1:nsub %in% bootsamp)){
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
