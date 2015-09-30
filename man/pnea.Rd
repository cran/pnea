\name{pnea}
\alias{pnea}
\alias{pvalue}
\alias{pneac}
\alias{netfromadj}
\title{Performs PNEA for lists of AGSs and FGSs}
\description{Compute parametric network enrichment analysis from altered gene sets to functional gene sets.}
\usage{
pnea(agslist, fgslist, network, nodes, alpha=NULL, agsnames, fgsnames)
}
\arguments{
  \item{agslist}{A list of altered gene sets (AGSs). Each element of the list is an AGS and must be represented as a vector of genes}
  \item{fgslist}{A list of functional gene sets (FGSs). Each element of the list is a FGS and must be represented as a vector of genes. This argument is optional: if provided, every set of \code{fgslist} is compared with every set of \code{agslist}; if \code{NULL}, the function compares sets in \code{agslist} between themselves (assuming that \code{network} is undirected)}
  \item{network}{It can be one of the following objects: an adjacency matrix (see 'Example 1'), an \code{igraph} object (see 'Example 2') or a two-column matrix where elements in the first column are parent nodes and in the second child nodes (see 'Example 3')}
  \item{nodes}{Vector containing the names of all nodes that form the network}
  \item{alpha}{Significance level of the test (optional). If specified, a column with the conclusion of the test is added to the output}
  \item{agsnames}{A vector of names for the AGSs (optional). It has to be provided only if the elements of agslist are not named}
  \item{fgsnames}{A vector of names for the FGSs (optional). It has to be provided only if the elements of fgslist are not named}
}
\value{A data frame with the following columns:
\item{AGS}{altered gene set}
\item{FGS}{functional gene set}
\item{naf}{observed number of links from AGS to FGS}
\item{expected_naf}{expected number of links from AGS to FGS (in absence of enrichment)}
\item{pvalue}{p-value of the test}
\item{conclusion}{conclusion of the test (only if \code{alpha} is specified): no enrichment, overenrichment or underenrichment}
}
\author{Mirko Signorelli}
\seealso{\code{\link{plot.pnea}}, \code{\link{print.pnea}}, \code{\link{summary.pnea}}}
\examples{
# Example 1: network given as adjacency matrix:
A = matrix(0, nrow=7, ncol=7)
A[1,c(2,3)]=1; A[2,c(5,7)]=1;A[3,c(1,4)]=1;A[4,c(2,5,7)]=1;A[6,c(2,5)]=1;A[7,4]=1
labels = letters[1:7]
set1 = c('a','e')
set2 = c('c','g')
set3 = c('d','f')
agslist = list('set 1' = set1, 'set 2' = set2)
fgslist = list('set 3' = set3)

test1 = pnea(agslist = agslist, fgslist = fgslist, network=A, nodes=labels, alpha=0.05)
print(test1)

agslist = list('set 1' = set1, 'set 2' = set2, 'set 3' = set3)
test2 = pnea(agslist = agslist, fgslist = NULL, network=A, nodes=labels, alpha=0.05)
print(test2)

# Example 2: network given as igraph object:
library(igraph)
igraphobject = erdos.renyi.game(15, 1/3)
set1 = c(1,4,7)
set2 = c(2,5)
set3 = c(3,9,14)
set4 = c(8,15)
agslist = list('set 1' = set1, 'set 2' = set2)
fgslist = list('set 3' = set3, 'set 4' = set4)

test2 = pnea(agslist, fgslist, network = igraphobject, nodes=seq(1,15), alpha=NULL)
print(test2)

# Example 3: network given as list of links:
networklist = matrix(nrow=12, ncol=2)
networklist[,1]=c('a','a','b','b','c','c','d','d','d','f','f','g')
networklist[,2]=c('b','c','e','g','a','d','b','e','g','b','e','d')
# REMARK: note the different ways used to specify directed and undirected edges. 
# Nodes 'a' and 'c' are connected by an undirected edge,
# whereas nodes 'b' and 'e' are linked by a directed arrow from 'b' to 'e'

labels = letters[1:7]
set1 = c('a','e')
set2 = c('c','g')
set3 = c('d','f')
agslist = list('set 1' = set1, 'set 2' = set2)
fgslist = list('set 3' = set3)

test3 = pnea(agslist, fgslist, network = networklist, nodes=labels, alpha=0.05)
print(test3)

# Example 4: ESR data
\dontrun{
data(yeast)
esr = list('ESR 1' = yeast$ags1, 'ESR 2' = yeast$ags2)
test = pnea(agslist = esr, fgslist = yeast$fgs, network = yeast$yeastnet,
nodes = yeast$ynetgenes, alpha = 0.01)

test1 = test[1:101,]
print(test1)
summary(test1)
plot(test1)
o1=test1[test1$conclusion=='Overenrichment',]
print(o1, nrows='ALL') #display overenrichments

test2 = test[102:202,]
print(test2)
summary(test2)
plot(test2, nbreaks=20)
o2=test2[test2$conclusion=='Overenrichment',]
print(o2, nrows='ALL') #display overenrichments
}}
\keyword{htest}
