
#library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R"); 

# Single-membership Stochastic Block Model: Reference Class

#input: ID.labels for a number of nodes, chosen to be the "dominant" ones.
#output: the permutation to switch labels to a simpler convention.

SBMcid <-
  setRefClass(
    "SBMcid",
    fields = list(
      n.groups="numeric",
      b.vector="numeric",

      b.vector.m="numeric",
      b.vector.v="numeric",
      
      
      membership="integer",
      #mult.factor="numeric",
      #mult.factor.m="numeric",
      #mult.factor.v="numeric",
      membership.a="matrix",

      shift="numeric",
      restrict.and.shift="logical",
      group.pairs="matrix",
      
      #inherited from main. Must fix later, but OK for now.
      n.nodes="numeric",
      outcome="numeric",
      edge.list="matrix",
      residual.variance="numeric",
      sr.rows="list"    #,
      ),
    
    methods=list(
      
      initialize = function (

        n.groups=1,
        
        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        sr.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),
       
        b.vector=rep(0, n.groups*(n.groups+1)/2),
        b.vector.m=rep(0, n.groups*(n.groups+1)/2),
        b.vector.v=rep(10000, n.groups*(n.groups+1)/2),
        
        membership=sample(n.groups, n.nodes, replace=TRUE),
        
        #mult.factor=1,
        #mult.factor.m=0,
        #mult.factor.v=1000000,

        membership.a=matrix(1, nrow=n.nodes, ncol=n.groups),
        shift=0,

        restrict.and.shift=TRUE,
        generate=FALSE
        
        ) {
        
        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$sr.rows <<- sr.rows
        
        .self$n.groups <<- n.groups
        
        .self$b.vector <<- b.vector
        .self$b.vector.m <<- b.vector.m
        .self$b.vector.v <<- b.vector.v
        .self$membership <<- membership

        #.self$mult.factor <<- mult.factor
        #.self$mult.factor.m <<- mult.factor.m
        #.self$mult.factor.v <<- mult.factor.v

        .self$membership.a <<- membership.a
        .self$residual.variance <<- residual.variance
        .self$restrict.and.shift <<- restrict.and.shift

        .self$group.pairs <<- makeEdgeListSelfies(n.groups)
                
        .self$shift <<- shift
        if (generate) .self$generate() else .self$outcome <<- outcome

        #center.me()
      },
      center.me = function () if (restrict.and.shift) {

        shift <<- mean(b.vector)
        b.vector <<- b.vector - shift
        
        #if (n.groups > 1) {
        #  sdbv <- sd(b.vector)
        #  b.vector <<- b.vector/sdbv
        #  mult.factor <<- mult.factor*sdbv
        #}
        
        #intercept <<- intercept + shift.t
      },
       
      reinitialize = function (n.nodes=NULL,
        edge.list=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          sr.rows <<- row.list.maker(edge.list)
        }
        if (length(membership) != n.nodes) {
          message ("Reinitializing SBM Memberships")
          membership <<- sample(n.groups, n.nodes, replace=TRUE)
          membership.a <<- matrix(1, nrow=n.nodes, ncol=n.groups)
        }
        
      },

      pieces = function (include.name=FALSE) {
        out <- list (b.vector=b.vector, membership=membership)#, mult.factor=mult.factor)
        class(out) <- "SBMout"
        #if (include.name) out <- c("SBM", out)
        out
      },

      show = function () {
        message("b.vector:"); print(b.vector)
        message("membership:"); print(membership)
        #message("mult.factor:"); print(mult.factor)
      },
      plot = function (memb=membership, block=symBlock(b.vector), ...) {
        single.membership.plot (memb, block, ...)
      },
      plot.network = function (color=outcome, ...) {
        netplot (edge.list, color, ...)
      },



      value = function () {
        sbm.matrix <- symBlock(b.vector)
        #mult.factor*
          sbm.matrix[membership[edge.list[,1]] +
                     dim(sbm.matrix)[1]*(membership[edge.list[,2]]-1)]
      },
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {   #slightly slower.
        sbm.matrix <- symBlock(parameters[[1]])
        #parameters[[3]]*
        sbm.matrix[parameters[[2]][edge.list[edges,1]] +
                   dim(sbm.matrix)[1]*(parameters[[2]][edge.list[edges,2]]-1)]
      },


      
      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},

      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },


      random.start = function () {
        membership <<- sample(n.groups, n.nodes, replace=TRUE)
        b.vector <<- rnorm(n.groups*(n.groups+1)/2, 0, 0.5)
        #if (restrict.and.shift) mult.factor <<- rnorm(1, mult.factor.m, sqrt(mult.factor.v))
      },
      
     # draw.mult.factor = function () {

     #   b.matrix <- symBlock(b.vector)
     #   X.term <- b.matrix[membership[edge.list[,1]] + (membership[edge.list[,2]]-1)*n.groups]
        
        #if (verbose>1) {print(X.term); print(outcome)}
          
     #   var.comp <- 1/(sum(X.term^2)/residual.variance + 1/mult.factor.v)
     #   mean.comp <- var.comp*(sum(X.term*outcome)/residual.variance + mult.factor.m/mult.factor.v)
     #   mult.factor <<- rnorm (1, mean.comp, sqrt(var.comp))

     # },

      rotate = function () {
        rotation <- SBM.ID.rotation(membership, n.groups)
        membership <<- rotation[membership]
        b.vector <<- SBM.rotate.bvector(b.vector, rotation)
      },
      
      draw = function (verbose=0) {

        if (length(outcome) != nrow(edge.list)) stop ("SBM: outcome and edge.list have different lengths.")

        b.matrix <- symBlock(b.vector)
        b.memb <- membership
        #b.factor <- mult.factor
        
        if (verbose>1) print(b.memb)

        # draw memberships.
        for (ii in sample(1:n.nodes)) {
          log.pp.vec <- sapply(1:n.groups, function(gg) {
            b.memb[ii] <- gg
            piece <- b.matrix[b.memb[edge.list[sr.rows[[ii]],1]] +
                              (b.memb[edge.list[sr.rows[[ii]],2]]-1)*n.groups]
            sum(dnorm(outcome[sr.rows[[ii]]], piece, sqrt(residual.variance), log=TRUE)) +
              log(membership.a[ii,gg])
          })
          log.pp.vec <- log.pp.vec - max(log.pp.vec)
          b.memb[ii] <- sample (1:n.groups, 1, prob=exp(log.pp.vec))
        }
        if (verbose>1) print(b.memb)
        membership <<- b.memb

        

        #draw block probs.
        #bits <- makeEdgeListSelfies(n.groups)
        membership.pairs <- cbind(b.memb[edge.list[,1]],
                                  b.memb[edge.list[,2]])
        
        b.vector <<- sapply(1:length(b.vector), function(bb) {
          picks <- unique(c(which(membership.pairs[,1]==group.pairs[bb,1] &
                                  membership.pairs[,2]==group.pairs[bb,2]),
                            which(membership.pairs[,1]==group.pairs[bb,2] &
                                  membership.pairs[,2]==group.pairs[bb,1])))
          if (length(picks) > 0) {
            var.b <- 1/(length(picks)/residual.variance + 1/b.vector.v[bb])
            mean.b <- var.b*(sum(outcome[picks])/residual.variance + b.vector.m[bb]/b.vector.v[bb])
            output <- rnorm(1, mean.b, sqrt(var.b))
          } else output <- rnorm(1, 0, 0.5)
          output
        })


        if (restrict.and.shift) {
          #message("Exterminate")
          center.me()
          #draw.mult.factor()
        }

        rotate()
        
      },

      gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1, make.random.start=FALSE) {
        out <- list()
        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
            if (report.interval > 0) if (index %% report.interval == 0) message("SBM ",index)
          }
        }
        return(out)
      },
      
      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.summary = function (gibbs.out) {
        membs1 <- {
          d1 <- sapply(gibbs.out, function(gg) number.to.vector(gg$membership))
          matrix(apply(d1, 1, mean),  ncol=n.nodes)
        }
        bvec <- apply(sapply(gibbs.out, function(gg) gg$b.vector), 1, mean)
        return(list(membership=membs1, b.vector=bvec))

      },
      
      gibbs.plot = function (gibbs.out, ...) {
        get.sum <- gibbs.summary(gibbs.out)
        block.membership.plot (get.sum[[1]], symBlock(get.sum[[2]]), main = "SBM Summary from Gibbs Sampler", ...)
      }


      )
    )

