\name{netfromadj}
\alias{netfromadj}
\title{
Converts an adjacency matrix into a matrix that can be processed by \code{pnea}}
\description{
Converts an adjacency matrix into a matrix that can be processed by \code{pnea}}
\usage{
netfromadj(A, nodes, nettype)
}
\arguments{
  \item{A}{The adjacency matrix representing the network}
  \item{nodes}{A vector with the names of nodes}
  \item{nettype}{Either \code{'directed'} or \code{'undirected'}}
}
\value{A two-column matrix, where every row represents and edge. For directed networks, elements in the first column must be parent nodes, and elements in the second column child nodes.}
\author{Mirko Signorelli}
\seealso{\code{\link{pnea}}}
