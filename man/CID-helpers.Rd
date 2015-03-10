\name{CIDnetworks-helpers}
\alias{l.diag}
\alias{u.diag}
\alias{ordinal.maker}
\alias{unwrap.CID.Gibbs}
\alias{mat.cov.to.edge.list.cov}
\alias{CID.metric}

\title{Functions to aid in the use and construction of CIDnetworks objects}
\description{
  Functions to aid in the use and construction of CIDnetworks objects
}
\usage{
  l.diag (nn)
  u.diag (nn)
  ordinal.maker (vec, cuts=quantile(vec, c(0.25, 0.5, 0.75)))

  unwrap.CID.Gibbs (gibbs.out)
  mat.cov.to.edge.list.cov (Xmat, n.nodes = dim(Xmat)[1],
                            arc.list = make.arc.list(n.nodes))

}
\arguments{

  \item{nn}{The number of rows in the square matrix for which we wish to
    extract the lower or upper diagonal matrix.}

  \item{vec}{The elements to be divided into ordinal categories.}
  \item{cuts}{The cut points at which to divide vec into ordinal
  categories. Default values separate vec into quartiles.}

  \item{gibbs.out}{The list object of draws from the Gibbs sampler. This
    re-sorts the object into a matrix form for easier consumption.}
  \item{Xmat}{A three-dimensional array of covariates, with n.nodes rows
    and columns. Each slice is a different covariate.}
  \item{n.nodes}{Number of nodes in network}
  \item{arc.list}{List of potential edges in network.}
  

}
\details{
  These functions are included for the convenience of users of
CIDnetworks.
  l.diag and u.diag provide the indices of a matrix to extract the lower and upper diagonal elements.
  ordinal.maker will turn any numeric vector into a series of ordinal
	integers for easy use in a CIDnetworks outcome.
	Xmat converts a sociomatrix-style array of covariates into one that
	can easily be used by the COV() component.

}
\author{A.C. Thomas <act@acthomas.ca>}

