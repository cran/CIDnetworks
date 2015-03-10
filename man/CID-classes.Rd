\name{CIDnetworks-classes}
\alias{BETA}
\alias{BETAcid}
\alias{EdgeCOV}
\alias{SenderCOV}
\alias{ReceiverCOV}
\alias{SendRecCOV}
\alias{IdenticalCOV}
\alias{COVARIATEcid}
\alias{HBM}
\alias{HBMcid}
\alias{LSM}
\alias{LSMcid}
\alias{LVM}
\alias{LVMcid}
\alias{MMSBM}
\alias{MMSBMcid}
\alias{SBM}
\alias{SBMcid}
\alias{SR}
\alias{SRcid}

\title{Classes of CIDnetwork subcomponents}
\description{
  Each of the reference classes used as components in the CIDnetworks methodology.
}
\usage{

  BETA(...)
  EdgeCOV (..., cov.type="Edge")
  SenderCOV (...)
  ReceiverCOV (...)
  SendRecCOV (...)
  IdenticalCOV (...)
  HBM (...)
  LSM (...)
  LVM (...)
  MMSBM (...)
  SBM (...)
  SR (...)

}
\arguments{

  \item{...}{Arguments passed to the respective classes.}
  \item{cov.type}{Specifies the type of covariate effect.  Edge
  specifies Edge-level covariate effects.}
    
}
\details{
  Each of these functions can serve as a subcomponent in the main
  CIDnetwork class object. Information on the nodes, edge list, and so
  forth will be passed automatically by any routine creating a CID
  object. Options are generally provided by default. Arguments for each class:

  BETA: required: (nothing). Parameters: intercept.sr.
  
  EdgeCOV: required: covariates (matrix). Parameters: Corresponding
  coefficient vector coef.cov.

  SenderCOV:  required: covariates (vector of length n.nodes).
              Parameters: Corresponding coefficient vector coef.cov.

  ReceiverCOV:  required: covariates (vector of length n.nodes).
                Parameters: Corresponding coefficient vector coef.cov.

  SendRecCOV:  required: covariates (vector of length n.nodes).
               Parameters: Corresponding coefficient vector coef.cov.

  IdenticalCOV:  required: covariates (vector of length n.nodes).
                 Parameters: Corresponding coefficient vector coef.cov.
	       

  HBM: required: n.groups (single value). Parameters: block.value,
  membership (for nodes to blocks), tree.parent (for blocks).

  LSM: required: dimension (single value). Parameters: latent.space.pos.

  LVM: required: dimension (single value). Parameters: latent.space.pos.

  MMSBM: required: n.groups (single value). Parameters: b.vector,
  membership.edge, membership.node.

  SBM: required: n.groups (single value). Parameters: b.vector,
  membership.

  SR: required: (nothing). Parameters: intercept.sr.

  
}
\value{
  Each expression yields a Reference Class object for the respective
  submodel. If generate=TRUE, it will produce an outcome value for that
  class depending on its specific properties.
}
\author{A.C. Thomas <act@acthomas.ca>}

