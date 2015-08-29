\name{summary.pneac}
\alias{summary.pneac}
\title{
Summary method of \code{pnea}
}
\description{
\code{summary} method for class "\code{pneac}"
}
\usage{
\method{summary}{pneac}(object, ...)
}
\arguments{
  \item{object}{An object of class "\code{pneac}"}
  \item{\dots}{Further arguments passed to or from other methods}
}
\value{
The \code{summary.pneac} function returns the following values:
\itemize{
\item the number of different altered gene sets (AGS) tested;
\item the number of different functional gene sets (FGS) tested;
\item the total number of comparisons;
\item the number of enrichments at 1\% and 5\% level;
\item the p-value of the Kolmogorov-Smirnov test to check if the distribution of p-values is uniform.
}
}
\author{Mirko Signorelli}
\seealso{
\code{\link{pnea}}, \code{\link{plot.pneac}}, \code{\link{summary.pneac}}
}
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
summary(test)
}
\keyword{methods}
