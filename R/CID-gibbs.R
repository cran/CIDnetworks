
#####################################################################################
#
# Gibbs Sampler collection for CID, given the collection of input terms.

#source("COV-reference.R"); source("SBM-reference.R"); source("LSM-reference.R"); source("SR-reference.R"); library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R"); 


# All subclasses available to CID.

LSM <- function(...) LSMcid$new(...)
LVM <- function(...) LVMcid$new(...)
SBM <- function(...) SBMcid$new(...)
SR <- function(...) SRcid$new(...)
COV <- function(...) COVcid$new(...)
MMSBM <- function(...) MMSBMcid$new(...)
HBM <- function(...) HBMcid$new(...)


unwrap.CID.Gibbs <- function (gibbs.out) list.output.to.matrices(gibbs.out)

#{  elements <- lapply (1:length(gibbs.out[[1]]), function(el) if (grepl("out", class(gibbs.out[[1]][[el]]))) {    r1 <- lapply(1:length(gibbs.out[[1]][[el]]), function (el2) sapply(gibbs.out, function(it) it[[el]][[el2]]))    names(r1) <- names(gibbs.out[[1]][[el]])    r1  } else sapply(gibbs.out, function(it) it[[el]]))  names(elements) <- names(gibbs.out[[1]])  elements}



CIDnetwork <-
  setRefClass (
    "CIDnetwork",
    fields = list(
      n.nodes="numeric",
      edge.list="matrix",
      sr.rows="list",
      outcome="numeric",
      node.names="character",
      
      class.outcome="character",
      ordinal.count="numeric",
      ordinal.cutoffs="numeric",

      int.outcome="numeric",
      
      intercept="numeric",
      intercept.m="numeric",
      intercept.v="numeric",
      
      residual.variance="numeric",
      residual.variance.ab="numeric",
      robit.augment="numeric",
      
      log.likelihood="numeric",
      
      components="list",
      comp.values="matrix"
      ),
    
    methods=list(
      initialize = function (
        
        edge.list,    #specs: this should be numbered from 1 to the max edge number.
        sociomatrix,
        
        sr.rows,
        n.nodes,      #this should come from the data.
        node.names=character(),
                
        intercept=0,
        intercept.m=0,
        intercept.v=1000000,
        
        residual.variance=1,
        residual.variance.ab=c(0.001, 0.001),
        
        outcome=numeric(0),
        generate=FALSE,

        class.outcome="ordinal",
        ordinal.count=2,
        ordinal.cutoffs=sort(rexp(ordinal.count-2, rate=50)),
        
        components=list()
        ) {

        if (!(class.outcome %in% c("binary","gaussian","ordinal"))) stop (paste("class.outcome",class.outcome,"is not supported. Must be one of", paste(c("binary","gaussian","ordinal"), collapse=",")))

        if (!missing(sociomatrix)) {
          if (!(class(sociomatrix) %in% c("array", "matrix"))) stop ("Sociomatrix must be a matrix.")
          if (nrow(sociomatrix) != ncol(sociomatrix)) stop ("Sociomatrix must be square.")
          new.nodes <- nrow(sociomatrix)
          .self$n.nodes <<- new.nodes
          .self$edge.list <<- make.edge.list (new.nodes)
          .self$outcome <<- sociomatrix[l.diag(new.nodes)]
        } else {
          if (missing(edge.list)) {
            if (missing(n.nodes)) stop ("Not detected: n.nodes, edge.list, sociomatrix") else {
              .self$n.nodes <<- n.nodes
              .self$edge.list <<- make.edge.list (n.nodes)
              .self$outcome <<- outcome
            }
          } else {
            if (missing(n.nodes)) .self$n.nodes <<- max(c(edge.list)) else .self$n.nodes <<- n.nodes
            .self$edge.list <<- edge.list
            .self$outcome <<- outcome
          }
        }

        if (!missing(sr.rows)) {
          .self$sr.rows <<- sr.rows
        } else {
          .self$sr.rows <<- row.list.maker(edge.list)
        }

        if (length(node.names) != .self$n.nodes) .self$node.names <<- as.character(1:.self$n.nodes) else .self$node.names <<- node.names
        

        if (class.outcome == "binary") {
          .self$class.outcome <<- "ordinal"
          .self$ordinal.count <<- 2
          .self$ordinal.cutoffs <<- numeric()
        } else {
          .self$class.outcome <<- class.outcome
          .self$ordinal.count <<- ordinal.count
          .self$ordinal.cutoffs <<- ordinal.cutoffs
        }

        .self$intercept <<- intercept
        .self$intercept.m <<- intercept.m
        .self$intercept.v <<- intercept.v
                
        if (class.outcome == "gaussian") {
          .self$residual.variance <<- residual.variance
        } else {
          .self$residual.variance <<- 1
        }
        .self$residual.variance.ab <<- residual.variance.ab
        

        #reinitialize components, just in case.
        if (length(components)>0) {
          if (class(components) != "list") components.t <- list(components) else components.t <- components
          for (kk in 1:length(components.t)) components.t[[kk]]$reinitialize (n.nodes, edge.list)
          .self$components <<- components.t
        } else .self$components <- list()
        
        #message("Component initialization complete.")
        
        if (generate) .self$generate() else if (length(outcome) == nrow(edge.list)) {
          .self$outcome <<- outcome
        } else stop (paste("In CIDnetwork$initialize: Outcome variable length: ",length(outcome),"; Edges: ", nrow(edge.list)))


        # Check for empty categories, currently disallowed.
        if (class.outcome == "ordinal") {
          counts <- sapply(0:(.self$ordinal.count-1), function(cc) sum(.self$outcome == cc))
          if (any(counts == 0)) stop ("The following ordinal categories have no outcomes: ", paste(which(counts == 0)-1, collapse=" "))
        }
        

        
        #message("Generation complete.")
        
        update.intermediate.outcome ()
        #message("CID Initialization complete.")
         
      },
      
      reinitialize = function (n.nodes=NULL,
        edge.list=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          sr.rows <<- row.list.maker(edge.list)
        }
        
        if (length(components)>0) for (kk in 1:length(components)) components[[kk]]$reinitialize (n.nodes, edge.list)
      },

      pieces = function (include.name=FALSE) {
        if (length(components)>0) {
          c(list(intercept=intercept,
                 residual.variance=residual.variance,
                 log.likelihood=log.likelihood,
                 ordinal.cutoffs=ordinal.cutoffs),
            lapply(components, function(cc) cc$pieces(include.name)))
        } else {
          list(intercept=intercept,
               residual.variance=residual.variance,
               log.likelihood=log.likelihood,
               ordinal.cutoffs=ordinal.cutoffs)
        }
      },

      show = function () {
        message("CIDnetwork object properties:")
        
        message(paste("class.outcome:", class.outcome))
        if (class.outcome == "ordinal") message(paste("Ordinal groups:", ordinal.count))
        message(paste("Nodes:", n.nodes))
        message(paste("Edges:", nrow(edge.list)))
        
        
        message("Intercept: ", intercept)
        message("Variance: ", residual.variance)
        if (length(components)>0) for (kk in 1:length(components)) {
          message(class(components[[kk]])); components[[kk]]$show()
        }

      },
      plot = function (coefs=coef.cov, names=1:length(coefs), sd=NULL, interval=NULL, ...) {
        if (length(components) > 0) for (cc in 1:length(components)) components[[cc]]$plot()
      },
      plot.network = function (color=outcome, ...) {
        netplot (edge.list, color, ...)
      },

      value = function (redo=FALSE) {
        if (redo) update.comp.values()
        rowSums(comp.values) + intercept
      },
      value.ext = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        if (length(components)>0) comp.values.here <- sapply (components, function(cc) cc$value()) else comp.values.here <- matrix(0, nrow=nrow(edge.list))
        (rowSums(comp.values.here) + parameters$intercept)[edges]
      },

      summary = function () show(),
      

      generate = function () {
        int.outcome <<- rnorm(nrow(edge.list),
                              value(redo=TRUE),
                              sqrt(residual.variance))
        if (class.outcome == "ordinal") {
          temp.out <- 0*int.outcome
          ordinal.steps <- c(0, ordinal.cutoffs)
          for (ii in 1:length(ordinal.steps)) temp.out <- temp.out + 1*(int.outcome > ordinal.steps[ii])
          outcome <<- temp.out
        }
        if (class.outcome == "gaussian") outcome <<- int.outcome        
      },
      
      
      rem.values = function(kk) {if (kk>0) value() - comp.values[,kk] else value() - intercept},
      

      

      update.intermediate.outcome = function () {
        if (class.outcome == "ordinal") {
          #first, update the values themselves.
          value.hold <- value(redo=TRUE)
          io.temp <- int.outcome

          breaker.lower <- c(-Inf, 0, ordinal.cutoffs)
          breaker.upper <- c(0, ordinal.cutoffs, Inf)
          
          #are any going to be trouble?
          p.gen <- rep(0.5, length(outcome))
          for (ii in 1:ordinal.count) {
            p.gen[outcome == ii-1] <-
              pnorm(breaker.upper[ii], value.hold[outcome == ii-1], 1) - pnorm(breaker.lower[ii], value.hold[outcome == ii-1], 1)
          }

          for (ii in 1:ordinal.count) {
            io.temp[outcome == ii-1 & p.gen > 1e-10] <-
              rtnorm(sum(outcome == ii-1 & p.gen > 1e-10),
                     value.hold[outcome == ii-1 & p.gen > 1e-10], 1,
                     lower=breaker.lower[ii], upper=breaker.upper[ii])
            io.temp[outcome == ii-1 & p.gen <= 1e-10] <- (breaker.lower[ii]+breaker.upper[ii])/2
            if (ii == 1) io.temp[outcome == ii-1 & p.gen <= 1e-10] <- breaker.upper[ii]
            if (ii == ordinal.count) io.temp[outcome == ii-1 & p.gen <= 1e-10] <- breaker.lower[ii]
          }
            
          int.outcome <<- io.temp

          #now, change the cutoff values, which lie between the Z values for each one. Assume a flat prior for now.
          if (length(ordinal.cutoffs)>0) for (kk in 1:length(ordinal.cutoffs)) {
            effective.range <- c(max(int.outcome[outcome <= kk]), min(int.outcome[outcome >= kk+1]))
            if (is.na(effective.range[1])) effective.range[1] <- 0
            if (is.na(effective.range[2])) effective.range[1] <- 10000            
            ordinal.cutoffs[kk] <<- runif(1, effective.range[1], effective.range[2])
          }
          
        }
        if (class.outcome == "gaussian") {int.outcome <<- outcome}
      },
      
      update.comp.values = function () {
        if (length(components)>0) comp.values <<- sapply (components, function(cc) cc$value()) else comp.values <<- matrix(0, nrow=nrow(edge.list))
      },

      draw.intercept = function (verbose=FALSE) {
        
        outcomeresid <- int.outcome - rem.values(0);
        
        varpiece <- solve(nrow(edge.list)/residual.variance + 1/intercept.v)
        meanpiece <- varpiece*(sum(outcomeresid)/residual.variance + intercept.m/intercept.v)
        if (verbose) message ("Intercept ",meanpiece," ",sqrt(varpiece))
        intercept <<- rnorm(1, meanpiece, sqrt(varpiece))

      },

      log.likelihood.by.value = function (value.this=value(), sumup=TRUE) {
        output <- NULL
        if (class.outcome == "gaussian") {
          outcomeresid <- int.outcome - value.this
          output <- dnorm(outcomeresid, 0, sqrt(residual.variance), log=TRUE)
        }
        
        if (class.outcome == "ordinal") {

          breaker.lower <- c(-Inf, 0, ordinal.cutoffs)
          breaker.upper <- c(0, ordinal.cutoffs, Inf)
          
          output <- 0*int.outcome
          for (ii in 1:(length(breaker.lower)-1))
            output[outcome == ii-1] <- log(
                     pnorm(rep(breaker.upper[ii], sum(outcome == ii-1)),
                           value.this[outcome == ii-1],
                           sqrt(residual.variance)) -
                     pnorm(rep(breaker.lower[ii], sum(outcome == ii-1)),
                           value.this[outcome == ii-1],
                           sqrt(residual.variance)))
                    
          #outs <- value.this
          #outs[outcome==0] <- -outs[outcome==0]
          #output <- pnorm(outs, 0, sqrt(residual.variance), log=TRUE)
        }
        if (sumup) output <- sum(output)
        return(output)
      },
      
      update.log.likelihood = function () {
        log.likelihood <<- log.likelihood.by.value ()        
      },
      
      draw.variance = function (verbose=FALSE) {
        outcomeresid <- int.outcome - value();

        if (verbose) message ("Variance: ",nrow(edge.list), " ", sum(outcomeresid^2))
        residual.variance <<-
          1/rgamma(1,
                   residual.variance.ab[1] + nrow(edge.list)/2,
                   residual.variance.ab[2] + sum(outcomeresid^2)/2)
        
      },
      
      draw = function (verbose=FALSE) {
        
        if (class.outcome != "gaussian") update.intermediate.outcome()
  
        if (length(components)>0) for (kk in 1:length(components)) {
          update.comp.values()
          components[[kk]]$outcome <<- .self$int.outcome - rem.values(kk)
          components[[kk]]$residual.variance <<- residual.variance
          components[[kk]]$draw()

          if (exists("shift", components[[kk]])) {
            intercept <<- intercept + components[[kk]]$shift
            components[[kk]]$shift <<- 0
          }
        }

        #variance and intercept.
        update.comp.values()
        draw.intercept(verbose)
        
        if (class.outcome == "gaussian") {
          update.comp.values()
          draw.variance(verbose)
        }
        
        update.log.likelihood()
        
      },
      
      random.start = function () {
        intercept <<- rnorm (1, 0, 1)
#        draw.intercept()
        if (length(components)>0) for (kk in 1:length(components)) components[[kk]]$random.start()
        if (class.outcome == "gaussian") draw.variance()
        update.log.likelihood()

      },
      
      gibbs.full = function (report.interval=100, draws=100, burnin=0, thin=1,
        make.random.start=TRUE) {
        
        out <- list()
        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- pieces(include.name=TRUE)
            if (report.interval > 0) if (index %% report.interval == 0) message("CID ",index)
          } else if (round(index)==index) {
            if (report.interval > 0) if (index %% report.interval == 0) message("CID burnin ",index)
          }
        }
        return(out)
      },

      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.switcheroo = function (gibbs.out) {
        out <- lapply(1:length(gibbs.out[[1]]), function(el)
                      lapply(1:length(gibbs.out), function(el2) gibbs.out[[el2]][[el]]))
        out
      },
      
      gibbs.summary = function (gibbs.out) {
        switched <- gibbs.switcheroo (gibbs.out)
        s.sum <- function (int1) c(min=min(int1), max=max(int1), mean=mean(int1), sd=sd(int1), quantile(int1, c(0.025, 0.975)))
        out <- list()

        out$intercept <- s.sum(unlist(switched[[1]]))
        out$residual.variance <- s.sum(unlist(switched[[2]]))
        out$log.likelihood <- s.sum(unlist(switched[[3]]))
        out$ordinal.cutoffs <- {
          "To be added"
        }
        
        if (length(components) > 0) for (cc in 1:length(components)) {
          out[[cc+4]] <- components[[cc]]$gibbs.summary(switched[[cc+4]])
          names(out)[cc+4] <- class(components[[cc]])
        }
          

        return(out)
      },

      gibbs.plot = function (gibbs.out, DIC=NULL, which.plots=1:(length(components)+4)) {
        switched <- gibbs.switcheroo (gibbs.out)
        
        if (1 %in% which.plots) plot.default (unlist(switched[[1]]), main="Grand Intercept")
        if (2 %in% which.plots & class.outcome=="gaussian") plot.default (unlist(switched[[2]]), main="Residual Variance")
        main.label <- "Log-likelihood"; if (!is.null(DIC)) main.label <- paste0(main.label, ": DIC = ",signif(DIC, 5))
        if (3 %in% which.plots)
          plot.default (unlist(switched[[3]]), main=main.label)
        if (4 %in% which.plots & class.outcome=="ordinal" & ordinal.count>2) {
          draws <- length(unlist(switched[[1]]))
          xx <- sort(rep(1:draws, ordinal.count-2))
          plot.default (c(1,xx), c(0,unlist(switched[[4]])), col=c(0, rep(1:(ordinal.count-2), draws)),
                        main="Ordinal Cutoff Values")
          abline(h=0, col=8)
        }
        
        
        if (length(components) > 0) for (cc in 1:length(components)) if (4+cc %in% which.plots) components[[cc]]$gibbs.plot(switched[[cc+4]])
        
      },
      
      DIC = function (gibbs.out) {
        all.values <- gibbs.value(gibbs.out)
        deviance.of.average <- -2*log.likelihood.by.value (apply(all.values, 1, mean))
        average.deviance <- mean(-2*apply(all.values, 2, log.likelihood.by.value))
        return(2*average.deviance - deviance.of.average)
      },

      marginal.loglikelihood = function (gibbs.out) {
        all.values <- gibbs.value(gibbs.out)
        model.log.likelihoods <- apply(all.values, 2, log.likelihood.by.value)
        1/mean(1/model.log.likelihoods)
      },

      pseudo.CV.loglikelihood = function (gibbs.out) {
        all.values <- gibbs.value(gibbs.out)
        model.log.likelihoods <- apply(all.values, 2, log.likelihood.by.value, sumup=FALSE)
        each.like <- log(1/apply(1/exp(model.log.likelihoods), 1, mean))
        sum(each.like)
      }

    )
    )

CID <- function (...) CIDnetwork$new(...)
CID.generate <- function (...) CIDnetwork$new(...)

CID.Gibbs <- function (edge.list,
                       outcome,
                       sociomatrix,

                       CID.object,
                       
                       #n.nodes=max(edge.list),
                       components=list(),
                       class.outcome=NULL,
                       fill.in.missing.edges=missing(outcome),
                         
                       ...) {
  #edge.list=n0$edge.list; outcome=n0$outcome; components=list(); n.nodes=max(edge.list)
  #edge.list=dolphins; components=list(LSM(2)); class.outcome=NULL

  if (missing(CID.object)) {
  
    if (missing(sociomatrix)) {
      if (missing(edge.list)) stop ("Neither an edge list not sociomatrix was provided.")

      edge.list <- as.matrix(edge.list)
      
      node.names <- unique(c(edge.list))
      n.nodes <- length(node.names)
      numbered.edge.list <- matrix (match(c(edge.list), node.names), ncol=2)
        
      if (missing(outcome) | fill.in.missing.edges) {
        if (missing(outcome)) message("Assuming that this is a complete network with specified edges as binary ties.") else message("Filling in unspecified edges as zeroes.")

        new.edge.list <- make.edge.list (n.nodes)
        rowmatch <- sapply (1:nrow(edge.list),
                            function(rr) min(which((numbered.edge.list[rr,1] == new.edge.list[,1] &
                                                    numbered.edge.list[rr,2] == new.edge.list[,2]) |
                                                   (numbered.edge.list[rr,2] == new.edge.list[,1] &
                                                    numbered.edge.list[rr,1] == new.edge.list[,2]))))

          
        temp.outcome <- rep(0, nrow(new.edge.list))
        if (missing(outcome)) {
          temp.outcome[rowmatch] <- 1
        } else {
          temp.outcome[rowmatch] <- outcome
        }
        
        outcome <- temp.outcome
        edge.list <- new.edge.list

      } else {
        edge.list <- numbered.edge.list
      }
    } else {
      n.nodes <- nrow(sociomatrix)
      edge.list <- make.edge.list (n.nodes)
      outcome <- sociomatrix[l.diag(n.nodes)]
      if (is.null(colnames(sociomatrix))) node.names <- 1:n.nodes else node.names <- colnames(sociomatrix)
    }

    ordinal.count <- 2
    if (is.null(class.outcome)) {
      if (all(round(outcome)==outcome)) {
        class.outcome <- "ordinal"
      } else {
        class.outcome <- "gaussian"
      }
    }
    if (class.outcome=="ordinal") {
      ordinal.count <- max(outcome)+1
      message ("Fitting: ordinal outcome with ",ordinal.count," states.")
    } else {
      message ("Fitting: gaussian outcome.")
    }
  
    CID.object <- CIDnetwork$new(edge.list, n.nodes=n.nodes, components=components, outcome=outcome,
                                 class.outcome=class.outcome, ordinal.count=ordinal.count,
                                 node.names=as.character(node.names))
  } 
    
  res <- CID.object$gibbs.full(...)
  DIC <- CID.object$DIC(res)
  #marg.ll <- CID.object$marginal.loglikelihood(res)
  #pseudo.CV.ll <- CID.object$pseudo.CV.loglikelihood(res)

  output <- list(results=res,  #unwrap.CID.Gibbs(
                 CID.object=CID.object,
                 DIC=DIC)
  class(output) <- "CID.Gibbs"
  
  return(output) #,
   #           marg.ll=marg.ll,
   #           pseudo.CV.ll=pseudo.CV.ll))
  
}


print.CID.Gibbs <- function (x, ...) {
  x$CID.object$show(...)
}

summary.CID.Gibbs <- function (object, ...) {
  object$CID.object$gibbs.summary (object$results, ...)
}

plot.CID.Gibbs <- function (x, ...) {
  x$CID.object$gibbs.plot (x$results, x$DIC, ...)
}

    
network.plot <- function (x, fitted.values=FALSE, ...) {

  if (class(x) == "CID.Gibbs") {

    if (fitted.values) values <- apply(x$CID.object$gibbs.value(x$results), 1, mean) else values <- x$CID.object$outcome
    x$CID.object$plot.network(values, ...)

  }
  
}

