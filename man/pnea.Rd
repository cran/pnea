\name{pnea}
\alias{pnea}
\alias{pneatest}
\alias{pneac}
\alias{netfromadj}
\title{Performs PNEA for lists of AGSs and FGSs}
\description{Given k altered gene sets and k functional gene sets, this function computes kq parametric network enrichment analyses (one for each couple of AGS and FGS).}
\usage{
pnea(agslist, fgslist, network, nodes, alpha=NULL, agsnames, fgsnames)
}
\arguments{
  \item{agslist}{A list of altered gene sets (AGSs). Each element of the list is an AGS and must be represented as a vector of genes}
  \item{fgslist}{A list of functional gene sets (FGSs). Each element of the list is a FGS and must be represented as a vector of genes}
  \item{network}{It can be one of the following objects: an adjacency matrix, an \code{igraph} object or a two-column matrix where elements in the first column are parent nodes and in the second child nodes.}
  \item{nodes}{Vector containing the names of all nodes in the network}
  \item{alpha}{Significance level of the test (optional). If specified, a column with the conclusion of the test is added to the output}
  \item{agsnames}{A vector of names for the AGSs (optional). It has to be provided only in case the elements of agslist are not named}
  \item{fgsnames}{A vector of names for the FGSs (optional). It has to be provided only in case the elements of fgslist are not named}
}
\value{A data frame with the following columns:
\item{AGS}{altered gene set}
\item{FGS}{functional gene set}
\item{naf}{number of links from AGS to FGS}
\item{pvalue}{p-value of the test}
\item{conclusion}{conclusion of the test (optional, only if \code{alpha} is specified): no enrichment, overenrichment, underenrichment}
}
\author{Mirko Signorelli}
\seealso{\code{\link{plot.pneac}}, \code{\link{print.pneac}}, \code{\link{summary.pneac}}}
\examples{
# Ex. 1: network given as adjacency matrix:
A = matrix(0, nrow=7, ncol=7)
A[1,c(2,3)]=1; A[2,c(5,7)]=1;A[3,c(1,4)]=1;A[4,c(2,5,7)]=1;A[6,c(2,5)]=1;A[7,4]=1
labels = letters[1:7]
set1 = c('a','e')
set2 = c('c','g')
set3 = c('d','f')
agslist = list('set 1' = set1, 'set 2' = set2)
fgslist = list('set 3' = set3)

test1 = pnea(agslist, fgslist, network=A, nodes=labels, alpha=0.05)
print(test1)

# Ex. 2: network given as list of links:
networklist = matrix(nrow=12, ncol=2)
networklist[,1]=c('a','a','b','b','c','c','d','d','d','f','f','g')
networklist[,2]=c('b','c','e','g','a','d','b','e','g','b','e','d')

test2 = pnea(agslist, fgslist, network = networklist, nodes=labels, alpha=0.05)
print(test2)

# Ex. 3: network given as igraph object:
library(igraph)
igraphobject = erdos.renyi.game(15, 1/3)
set1 = c(1,4,7)
set2 = c(2,5)
set3 = c(3,9,14)
set4 = c(8,15)
agslist = list('set 1' = set1, 'set 2' = set2)
fgslist = list('set 3' = set3, 'set 4' = set4)

test3 = pnea(agslist, fgslist, network = igraphobject, nodes=seq(1,15), alpha=NULL)
print(test3)
}
\keyword{htest}
