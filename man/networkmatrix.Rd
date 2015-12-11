\name{networkmatrix}
\alias{networkmatrix}
\title{
Creates a network matrix for \code{pnea}
}
\description{
Internal function, creates a two-column network matrix that can be further processed by \code{\link{pnea}}.
}
\usage{
networkmatrix(network, nodes, nettype)
}
\arguments{
  \item{network}{
One of the following objects: an adjacency matrix (class \code{"matrix"}), a sparse adjacency matrix (class \code{"dgCMatrix"}) or an igraph graph (class \code{"igraph"})
}
  \item{nodes}{
Vector containing the (ordered) names of all nodes in the network
}
  \item{nettype}{
Either \code{'directed'} or \code{'undirected'}
}
}
\details{
This is an internal function, that is called within \code{\link{pnea}} to convert different types of network objects (see argument 'network' above) into a standard two-column network matrix, that can then be processed by \code{pnea}. 
}
\value{
A two-column matrix, where every row represents and edge. For directed networks, parent nodes must be in the first column, and child nodes in the second.
}
\author{Mirko Signorelli}
\seealso{\code{\link{pnea}}}
\examples{
# First case: adjacency matrix
n<-50
adjacency <- matrix(sample(0:1, n^2, replace=TRUE, prob=c(0.9,0.1)), ncol=n)
diag(adjacency) <- 0
lab = paste(rep('gene'),1:n)
head(networkmatrix(adjacency, lab, 'directed'))

# Second case: sparse adjacency matrix
library(Matrix)
sparse_adjacency<-Matrix(adjacency,sparse=TRUE)
head(networkmatrix(sparse_adjacency, lab, 'directed'))

# Third case: igraph object
library(igraph)
igraph_graph = erdos.renyi.game(15, 1/3)
lab = paste(rep('gene'),1:15)
head(networkmatrix(igraph_graph, lab, 'directed'))
}
\keyword{manip}