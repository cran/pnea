\name{summary.pnea}
\alias{summary.pnea}
\title{
Summary method of \code{pnea}
}
\description{
\code{summary} method for class "\code{pnea}".
}
\usage{
\method{summary}{pnea}(object, ...)
}
\arguments{
  \item{object}{An object of class "\code{pnea}"}
  \item{\dots}{Further arguments passed to or from other methods}
}
\value{
The \code{summary.pnea} function returns the following values:
\itemize{
\item the number of tests computed;
\item the number of enrichments at 1\% and 5\% level;
\item the p-value of the Kolmogorov-Smirnov test to check if the distribution of p-values is uniform.
}
}
\author{Mirko Signorelli}
\seealso{
\code{\link{pnea}}, \code{\link{plot.pnea}}, \code{\link{summary.pnea}}
}
\examples{
\dontrun{
data(yeast)
esr = list('ESR 1' = yeast$esr1, 'ESR 2' = yeast$esr2)
test = pnea(alist = esr, blist = yeast$goslimproc, network = yeast$yeastnet,
            nettype = 'undirected', nodes = yeast$ynetgenes, alpha = 0.01)

test1 = test[1:99,]
summary(test1)

test2 = test[100:198,]
summary(test2)
}}
\keyword{methods}
