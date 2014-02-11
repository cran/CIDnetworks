
#library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R"); 
# Standard Latent Space Model: Reference Class


LSMcid <-
  setRefClass(
    "LSMcid",
    fields = list(
      dimension="numeric",
      
      latent.space.pos="matrix",
      mult.factor="numeric",
      
      #mult.factor.m="numeric",
      #mult.factor.v="numeric",
      
      latent.space.pos.m="numeric",
      latent.space.pos.v="numeric",
      tune="numeric",

      #inverted.model="logical",
                     #inherited from main. Must fix later, but OK for now.

      node.names="character",
      n.nodes="numeric",
      outcome="numeric",
      edge.list="matrix",
      residual.variance="numeric",
      sr.rows="list"    #,
      
      ),
    
    methods=list(
      
      initialize = function (
        
        dimension=1,
        
        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        sr.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),
        
        latent.space.pos=matrix(rnorm(dimension*n.nodes), nrow=n.nodes),
        #mult.factor.m=0,
        #mult.factor.v=10000,
        
        latent.space.pos.m=0,
        latent.space.pos.v=10000,
        tune=0.1,

        inverted.model=FALSE,
#        mult.factor=-1 + 2*inverted,

        generate=FALSE
        
        ) {
        
        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$sr.rows <<- sr.rows
        .self$residual.variance <<- residual.variance
        .self$node.names <<- as.character(1:.self$n.nodes)
                
        .self$dimension <<- dimension
        .self$latent.space.pos <<- latent.space.pos
        
        #.self$mult.factor <<- mult.factor
        #.self$mult.factor.m <<- mult.factor.m
        #.self$mult.factor.v <<- mult.factor.v
        .self$latent.space.pos.m <<- latent.space.pos.m
        .self$latent.space.pos.v <<- latent.space.pos.v
        .self$mult.factor <<- -1 + 2*inverted.model
       
        .self$tune <<- tune

        #adjust.lsp()
        
        if (generate) .self$generate() else .self$outcome <<- outcome

      },

      #adjust.lsp = function (mult.up=TRUE) {
      #  mft <- mean(edge.list.distance(latent.space.pos, edge.list))
      #  mult.factor <<- mult.factor*mft
      #  latent.space.pos <<- latent.space.pos/mft          
      #},
      
      reinitialize = function (n.nodes=NULL,
        edge.list=NULL, node.names=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes  #.self$
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          sr.rows <<- row.list.maker(edge.list)
        }
        if (nrow(latent.space.pos) != .self$n.nodes) {
          message ("Reinitializing LSM Positions")
          latent.space.pos <<- matrix(rnorm(dimension*n.nodes), nrow=n.nodes)
       #   adjust.lsp()
        }
        if (!is.null(node.names)) {
          if (length(node.names) == .self$n.nodes) node.names <<- node.names
        } else node.names <<- as.character(1:.self$n.nodes)        
      },
      
      pieces = function (include.name=FALSE) {
        out <- list (latent.space.pos=latent.space.pos, mult.factor=mult.factor)
        class(out) <- "LSMout"
        out
      },
      
      show = function () {
        message("t(latent.space.pos):"); print(t(latent.space.pos))
        message("mult.factor:"); print(mult.factor)
      },
      plot = function (pos=latent.space.pos, ...) {
        latent.space.plot (pos, labels=node.names, ...)
      },
      plot.network = function (color=outcome, ...) {
        netplot (edge.list, color, node.labels=node.names, ...)
      },

      
      value = function () {mult.factor*edge.list.distance(latent.space.pos, edge.list)},
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {   #slightly slower.
        parameters[[2]]*edge.list.distance(parameters[[1]], rbind(edge.list[edges,]))
      },
      

      
      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},
      
      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },

      
      random.start = function () {
        latent.space.pos <<- matrix(rnorm(dimension*n.nodes), nrow=n.nodes)
        #mult.factor <<- rnorm(1, mult.factor.m, sqrt(mult.factor.v))
      },
      
      draw = function (verbose=0, mh.tune=tune) {
      #d1 <- LSMcid$new(); latent.space.pos <- d1$latent.space.pos; mult.factor <- d1$mult.factor; edge.list <- d1$edge.list; sr.rows <- d1$sr.rows; n.nodes <- d1$n.nodes; mult.factor.m=0; mult.factor.v=10000; tune=0.1

        lsdim <- dim(latent.space.pos)[2]

        latent.space.pos.hold <- latent.space.pos
        latent.space.pos.prop <- latent.space.pos
        #mult.factor.hold <- mult.factor
        
        for (dd in 1:n.nodes) {
                                        #proposal for latent space positions.
          latent.space.pos.prop[dd,] <- latent.space.pos.prop[dd,] +
            c(rmvnorm(1,
                      rep(0, lsdim),
                      diag(mh.tune^2, lsdim)))
          log.dens.prop <- log.likelihood(list(latent.space.pos.prop, mult.factor), sr.rows[[dd]])
          log.dens.orig <- log.likelihood(list(latent.space.pos.hold, mult.factor), sr.rows[[dd]])
          
          if (log.dens.prop - log.dens.orig > -rexp(1))
            latent.space.pos.hold[dd,] <- latent.space.pos.prop[dd,] else
          latent.space.pos.prop[dd,] <- latent.space.pos.hold[dd,]
          
          
           #correct scale factor. simple normal after adjustment!
        #  latent.space.pos.hold <- latent.space.pos.hold/mean(edge.list.distance(latent.space.pos.hold, edge.list))
        #  distances <- edge.list.distance(latent.space.pos, edge.list)
        #  cls.var <- 1/(distances %*% distances / residual.variance + 1/mult.factor.v)
        #  cls.mean <- (distances %*% outcome/residual.variance + mult.factor.m/mult.factor.v)*cls.var
        #  mult.factor.hold <- rnorm(1, cls.mean, sqrt(cls.var))
          
        }

  #Rotate back.
        latent.space.pos.hold <- postprocess.latent.positions(latent.space.pos.hold)
        rownames(latent.space.pos.hold) <- node.names
        
        latent.space.pos <<- latent.space.pos.hold
        #mult.factor <<- mult.factor.hold
        
  #return(list(latent.space.pos=latent.space.pos, coef.latent.space=coef.latent.space))
 
      },
      
      gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1, make.random.start=FALSE) {
        out <- list()
        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
            if (report.interval > 0) if (index %% report.interval == 0) message("LSM ",index)
          } else if (round(index)==index) {
            if (report.interval > 0) if (index %% report.interval == 0) message("LSM burnin ",index)
          }
        }
        return(out)
      },

      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.summary = function (gibbs.out) {
        lsp.all <- sapply(gibbs.out, function(gg) gg$latent.space.pos)
        output <- matrix(apply(lsp.all, 1, mean), nrow=n.nodes)
        rownames(output) <- node.names
        return(output)
      },
      
      gibbs.plot = function (gibbs.out, ...) {
        get.sum <- gibbs.summary(gibbs.out)
        plot (get.sum, main = "Mean Latent Space Positions from Gibbs Sampler", ...)
      }
      
      )
    )

  
  

