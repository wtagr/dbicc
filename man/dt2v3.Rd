\name{dt2v3}
\alias{dt2v3}
\docType{data}
\title{
Dataset: Functional Connectivity in the Human Brain
}
\description{
The resting-state fMRI data set, includes 25 participants (mean age 29.44 +/- 8.64, 10 males) scanned at New York University.
}
\usage{data("dt2v3")}
\format{
  The format is: List of 50 including 25 individuals * 2 scans

}
\details{
A Siemens Allegra 3.0-Tesla scanner was used to obtain three resting-state scans
for each participant drawn out the second and third scans in this dataset. Each scan consisted of 197 contiguous EPI functional volumes. The data were preprocessed using the FSL (http://www.fmrib.ox.ac.uk) and AFNI (http://afni.nimh.nih.gov/afni) software packages. The ROIs for our connectivity analysis are derived from the work of Gordon et al. (2016).
}
\source{
http://www.nitrc.org/projects/nyu_trt
}
\references{
\itemize{
\item P.T. Reiss, M. Xu, I. Cribben (2019). Generalized reliability based on distances. \url{https://arxiv.org/abs/1912.07137}.

\item Gordon, E. M., Laumann, T. O., Adeyemo, B., Huckins, J. F., Kelley, W. M., and Petersen,
S. E. (2016). Generation and evaluation of a cortical area parcellation from resting-state
correlations. Cerebral Cortex 26, 288-303.
}
}
\examples{
# load data
data(dt2v3)
# output the distance matrix for dm2icc
mdist<-function(datalist,method=c("f")){
  dmat<-matrix(0,length(datalist),length(datalist))
  for (i in 2:length(datalist)) for (j in 1:(i-1)){
    if(method=="r"){
      v1<-datalist[[i]][lower.tri(datalist[[i]], diag = FALSE)]
      v2<-datalist[[j]][lower.tri(datalist[[j]], diag = FALSE)]
      r<-cor(v1,v2)
      dmat[i,j]=dmat[j,i]=sqrt(2)*sqrt(1-r)
    }else if (method=='r/2'){
      v1<-datalist[[i]][lower.tri(datalist[[i]], diag = FALSE)]
      v2<-datalist[[j]][lower.tri(datalist[[j]], diag = FALSE)]
      r<-cor(v1,v2)
      dmat[i,j]=dmat[j,i]=sqrt(1-r)
    }else if(method=='sl1'){
      dmat[i,j]=dmat[j,i]=sum(abs(datalist[[i]]-datalist[[j]]))
    }else if(method=='sl4'){
      dmat[i,j]=dmat[j,i]=(sum((datalist[[i]]-datalist[[j]])^4))^(1/4)
    }else{
      dmat[i,j]=dmat[j,i]=norm(datalist[[i]]-datalist[[j]],type=method)
    }
  }
  return(dmat)
}

#distance matrix for correlation matrices
cor2v3<-lapply(dt2v3,cor)
distmat<-mdist(cor2v3)

#plot distance matrix
plotdmat(distmat,25,2)

#calculate the reliability
dm2icc(distmat,25,2)

}
\keyword{datasets}
