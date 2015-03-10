\name{CIDnetworks-master}
\alias{CID}
\alias{CID.generate}
\alias{CIDnetwork}

\title{Forward Simulation and Declaration of CID Network Models}
\description{
  Usage of the CID main class.
}
\usage{

  CID (input, outcome, n.nodes, node.names, intercept = 0,components,
  class.outcome="ordinal", fill.in.missing.edges=missing(outcome),
  generate=FALSE, verbose=2, ...)
  CID.generate (...)

}
\arguments{
  
  \item{input}{A matrix object containing information about the edges in a
    network.  If the input is a matrix with 2 columns the input is
    assumed to be a set of edges.  If input is a square matrix, it is
    assumed to be a sociomatrix.  If no input is provided, a new network
    is generated containing n.nodes.}
  \item{outcome}{A set of outcome values for each edge in the edgelist
    given as input.  If outcome is missing, the edgelist is assumed to
    represent observed edges in a binary network.  outcome requires no
    value is a sociomatrix is provided.}
  \item{n.nodes}{Number of nodes in generated networks}
  \item{node.names}{Names labeling each node in the network.}
  \item{intercept}{Intercept to be used when generating new networks}
  \item{components}{A list of components to use when generating new
    networks.  Possible components include SR, SBM, MMSBM, LSM, LVM HBM
    EdgeCOV, SenderCOV, ReceiverCOV, SendRecCov, IdenticalCOV}
  \item{class.outcome}{One of "ordinal" (default, values from 0 to
    higher integers), "binary" (ordinal in 0 and 1) or "gaussian"
    (unbounded continuous values). Class is auto-detected if NULL remains
    in place.}
  \item{fill.in.missing.edges}{If TRUE, the edge list will be augmented
    with zeroes for all unspecified but possible edges. By default, if an
    outcome is specified, these edges will not be added.}
  \item{generate}{When TRUE a new network with n.nodes nodes is
    generated as specified by the intercept and components passed to CID.}
  \item{verbose}{Level of output to be displayed while running.  A value
    of 0 will return little or no output.  A value of 1 will only notify
    of warnings of misuse.  A value of 2 will report progress of MCMC
    chain.  A value of 3 or higher will report debugging values.}

  \item{...}{Arguments passed to the master class.}
    
}
\details{
  CID generates a CIDNetwork object that can be passed to many of the
  plotting and examination methods in the CIDNetworks package.  Initial
  conversion of network data in standard edgelist or sociomatrix form
  allows easier manipulation within the CIDNetworks framework.

  CID can also be used to simulate networks from any of the models that
  CIDNetworks can fit.  CID.generate is a simple wrapper which calls CID
  with generate = TRUE.
}
\value{
  A CIDnetwork object.
}
\author{A.C. Thomas <act@acthomas.ca>}

