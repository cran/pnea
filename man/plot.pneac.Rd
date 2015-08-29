\name{plot.pneac}
\alias{plot.pneac}
\title{
Plot method of \code{pnea}
}
\description{
Histogram and p-p plot of p-values for class "\code{pneac}".
}
\usage{
\method{plot}{pneac}(x, nbreaks = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "\code{pneac}"}
  \item{nbreaks}{Number of breaks to be used for the histogram (default is 10)}
  \item{\dots}{Further arguments passed to or from other methods}
}
\value{
An histogram showing the distribution of p-values and a p-p plot comparing the distribution of p-values to the uniform distribution.}
\author{Mirko Signorelli}
\seealso{\code{\link{pnea}}, \code{\link{print.pneac}}, \code{\link{summary.pneac}}}
\examples{
A = matrix(0, nrow=7, ncol=7)
A[1,c(2,3)]=1; A[2,c(5,7)]=1;A[3,c(1,4)]=1;A[4,c(2,5,7)]=1;A[6,c(2,5)]=1;A[7,4]=1

labels = letters[1:7]
set1 = c('a','e')
set2 = c('c','g')
set3 = c('d','f')
agslist = list('set 1' = set1, 'set 2' = set2)
fgslist = list('set 3' = set3)

test = pnea(agslist, fgslist, network=A, nodes=labels, alpha=0.05)
plot(test)
}
\keyword{methods}