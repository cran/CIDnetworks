\name{CIDnetworks-gibbs}
\alias{CID.Gibbs}
\alias{network.plot}
\alias{sociogram.plot}
\alias{print.CID.Gibbs}
\alias{summary.CID.Gibbs}
\alias{print.summary.CID.Gibbs}
\alias{plot.CID.Gibbs}
\alias{likelihood.plot}
\alias{intercept.plot}
\alias{COV.plot}
\alias{LSM.plot}
\alias{SBM.plot}
\alias{MMSBM.plot}
\alias{SR.plot}
\alias{n.nodes}
\alias{edge.list}
\alias{outcome}
\alias{is.net.directed}
\alias{net.density}
\alias{node.names}
\alias{inDegree}
\alias{outDegree}
\alias{socio}
\alias{value.mat}
\alias{value.mat.mean}
\alias{switcheroo}



\title{CIDnetworks: Bayesian Inference}
\description{
  Usage of CID Gibbs Samplers.
}
\usage{

  CID.Gibbs (input,
             outcome,

             components,
             class.outcome = NULL,
             fill.in.missing.edges = missing(outcome),
             new.chain = FALSE,

             draws = 100,
             burnin = -1,
             thin = 10,
             report = 100,
             auto.converge = FALSE,
             extend.max=10,
             extend.count=100,
             verbose=2,
             ...)

  \method{print}{CID.Gibbs} (x, \dots)
  \method{summary}{CID.Gibbs} (object, \dots)
  \method{plot}{CID.Gibbs} (x, \dots)
  \method{print}{summary.CID.Gibbs} (x, \dots)

  likelihood.plot(x, ...)
  intercept.plot(x, mode = c("standard","trace"), ...)
  COV.plot(x, mode = c("standard","trace","scatterplot"), ...)
  LSM.plot(x, ...)
  SBM.plot(x, ...)
  MMSBM.plot(x, ...)
  SR.plot(x, ...)

  network.plot (x, fitted.values=FALSE, ...)
  sociogram.plot (x, component.color=0, vertexcolor, add.labels = TRUE, ...)

  n.nodes(object)
  edge.list(object)
  is.net.directed(object)
  net.density(object)
  outcome(object)
  node.names(object)
  inDegree(object)
  outDegree(object)
  socio(object)
  value.mat(CID.Gibbs.object, prob = TRUE)
  value.mat.mean(object, prob = TRUE)
  switcheroo(CID.Gibbs.object)  
 
  

}
\arguments{

  \item{input}{An object containing information about the edges in a
    network.  Must be one of the following classes:  Matrix, CIDnetwork, or
    CID.Gibbs.  If input is a square matrix, it is assumed to be a
    sociomatrix.  Otherwise a matrix with 2 columns and a number of rows
    equal to the number of edges is required.  Providing a CIDnetwork will
    use the associated edge.list.  Providing a CID.Gibbs object will
    continue the MCMC chain from the last draw.}

  \item{outcome}{If an edgelist is provided as input and outcome is
    missing, the edges provided are assumed to be the ones in a binary
    sociomatrix.  Otherwise, an outcome value must be specified for each
    edge in the edgelist, and any edges not provided are assumed to have
    no data.}
  
  \item{CID.Gibbs.object}{If desired, an existing CID.Gibbs output object can be loaded
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
  \item{new.chain}{If a CID.Gibbs object is provided, the default value
    of FALSE will return both the old and new MCMC chain combined.  A
    value of TRUE will drop the old chain completely.}
  \item{draws}{Number of draws to return.}
  \item{burnin}{Number of draws to burnin.  A negative value will
    automatically determine burnin amount.}
  \item{thin}{Amount of draws to thin the chain by.}
  \item{report}{Number of draws between reporting total draws so far.}
  \item{auto.converge}{When true, a Geweke convergence test on
    log-likleihood to detect convergence.}
  \item{extend.max}{Maximum number of times chain will be extended until
    it returns without converging.}
  \item{extend.count}{Number of draws to extend chain by if convergence
    test fails}
  \item{verbose}{Level of output to be displayed while running.  A value
    of 0 will return little or no output.  A value of 1 will only notify
    of warnings of misuse.  A value of 2 will report progress of MCMC
    chain.  A value of 3 or higher will report debugging values.}
  
  \item{...}{Further arguments to be passed to the Gibbs sampler
    routine or the plot routine. See details for more.}
  
  \item{x, object}{An object outputted from CID.Gibbs.}
  \item{mode}{Controls which diagnostic plot is made.}
  \item{fitted.values}{If TRUE, plots the fitted tie strength under the
    Gibbs sampler. If FALSE, plots the network outcomes as entered.}
  \item{component.color}{If non-zero, colors the nodes in the sociogram
    according to the output of the Gibbs sampler.}
  \item{vertexcolor}{User-passed vertex colors for sociogram.plot .}
  \item{add.labels}{When true, node labels are included on nodes.}
	\item{trace}{If selected, displays the Gibbs sampler trace plot for
	the intercept rather than a point and interval.}
	\item{prob}{In value.mat, converts the linear predictor value to the
	probability of a binary edge.}
     
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
