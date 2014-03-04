\name{CIDnetworks-gibbs}
\alias{CID.Gibbs}
\alias{network.plot}
\alias{sociogram.plot}
\alias{unwrap.CID.Gibbs}
\alias{print.CID.Gibbs}
\alias{plot.CID.Gibbs}
\alias{summary.CID.Gibbs}

\title{CIDnetworks: Bayesian Inference}
\description{
  Usage of CID Gibbs Samplers.
}
\usage{

  CID.Gibbs (edge.list,
             outcome,
             sociomatrix,

             CID.object,

             components=list(),
             class.outcome=NULL,
             fill.in.missing.edges=missing(outcome),
                         
             ...)
  unwrap.CID.Gibbs (gibbs.out)

  \method{print}{CID.Gibbs} (x, \dots)
  \method{summary}{CID.Gibbs} (object, \dots)
  \method{plot}{CID.Gibbs} (x, \dots)

  network.plot (x, fitted.values=FALSE, ...)
  sociogram.plot (x, component.color=0, ...)
 
  

}
\arguments{

  \item{edge.list}{A list of (potential) edges present in the network.}
  \item{outcome}{The outcome on each of these edges.}
  \item{sociomatrix}{If preferred, a sociomatrix of edges and
  values. The lower triangle will be taken in.}

  \item{CID.object}{If desired, an existing CID object can be loaded
  instead of a new network specification.}

  \item{components}{A list of sub-components, including (COV, HBM, LSM,
  LVM, MMSBM, SBM, SR).}
  \item{class.outcome}{One of "ordinal" (default, values from 0 to
  higher integers), "binary" (ordinal in 0 and 1) or "gaussian"
  (unbounded continuous values). Class is auto-detected if NULL remains
  in place.}
  \item{fill.in.missing.edges}{If TRUE, the edge list will be augmented
  with zeroes for all unspecified but possible edges. By default, if an
  outcome is specified, these edges will not be added.}
  
  \item{...}{Further arguments to be passed to the Gibbs sampler
  routine. See details for more.}

  \item{gibbs.out}{The list object of draws from the Gibbs sampler. This
    re-sorts the object into a matrix form for easier consumption.}

  \item{x, object}{An object outputted from CID.Gibbs.}
  \item{fitted.values}{If TRUE, plots the fitted tie strength under the
    Gibbs sampler. If FALSE, plots the network outcomes as entered.}
  \item{component.color}{If non-zero, colors the nodes in the sociogram
    according to the output of the Gibbs sampler.}
  
     
}
\details{
  This is the main routine for running a Gibbs sampler on any of the CID
  models. See the vignettes for more information.
}
\value{
  CID.Gibbs outputs a list containing a CID object, the results of the Gibbs sampler, and
  the Deviance Information Criterion estimate for the Gibbs.
}
\author{A.C. Thomas <act@acthomas.ca>}
