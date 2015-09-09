\name{plot.pnea}
\alias{plot.pnea}
\title{
Plot method of \code{pnea}
}
\description{
Histogram and p-p plot of p-values for objects of class "\code{pnea}".
}
\usage{
\method{plot}{pnea}(x, nbreaks = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "\code{pnea}"}
  \item{nbreaks}{Number of breaks to be used in the histogram (default is 10)}
  \item{\dots}{Further arguments passed to or from other methods}
}
\value{
An histogram showing the distribution of p-values and a p-p plot comparing the distribution of p-values to the uniform distribution.}
\author{Mirko Signorelli}
\seealso{\code{\link{pnea}}, \code{\link{print.pnea}}, \code{\link{summary.pnea}}}
\examples{
\dontrun{
data(yeast)
esr = list('ESR 1' = yeast$ags1, 'ESR 2' = yeast$ags2)
test = pnea(agslist = esr, fgslist = yeast$fgs, network = yeast$yeastnet,
nodes = yeast$ynetgenes, alpha = 0.01)

test1 = test[1:101,]
plot(test1)

test2 = test[102:202,]
plot(test2, nbreaks=20)
}}
\keyword{methods}