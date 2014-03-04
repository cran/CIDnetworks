
#library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R"); 

# Single-membership Stochastic Block Model: Reference Class

#input: ID.labels for a number of nodes, chosen to be the "dominant" ones.
#output: the permutation to switch labels to a simpler convention.

draw.MMSB.from.nodes <- function (edge.list, node.block) {
  block.id <- sapply(c(edge.list), function(jj)
                     sample(1:nrow(node.block), 1, prob=node.block[,jj]))
  matrix(block.id, ncol=2)
}

MMSBMcid <-
  setRefClass(
    "MMSBMcid",
    fields = list(
      n.groups="numeric",
      b.vector="numeric",
      b.vector.m="numeric",
      b.vector.v="numeric",

      membership.edge="matrix",   #looks like edge list: has the membership number of each participant.
      membership.node="matrix",   #each column is the Dirichlet distribution for each membership.
      membership.alpha0="numeric",  #prior strength for the Dirichlet -- alpha0*vector(1)
      
      shift="numeric",
      restrict.and.shift="logical",
      group.pairs="matrix",

      #single.membership="logical",
      
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

        n.groups=2,
        
        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        sr.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),
       
        b.vector=rep(0, n.groups*(n.groups+1)/2),
        b.vector.m=rep(0, n.groups*(n.groups+1)/2),
        b.vector.v=rep(10000, n.groups*(n.groups+1)/2),

        membership.alpha0=0.1,
        membership.node=rdirichlet.block (matrix(membership.alpha0, nrow=n.groups, ncol=n.nodes)),
        membership.edge=draw.MMSB.from.nodes(edge.list, membership.node),
        
        shift=0,
        #single.membership=FALSE,
        
        restrict.and.shift=FALSE,
        generate=FALSE
        
        ) {
        
        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$sr.rows <<- sr.rows
        .self$node.names <<- as.character(1:.self$n.nodes)
        
        .self$n.groups <<- n.groups
        
        .self$b.vector <<- b.vector
        .self$b.vector.m <<- b.vector.m
        .self$b.vector.v <<- b.vector.v
        .self$membership.edge <<- membership.edge
        .self$membership.node <<- membership.node
        .self$membership.alpha0 <<- membership.alpha0

        .self$residual.variance <<- residual.variance
        .self$restrict.and.shift <<- restrict.and.shift
        #.self$single.membership <<- FALSE

        .self$group.pairs <<- makeEdgeListSelfies(n.groups)
                
        .self$shift <<- shift
        if (generate) .self$generate() else .self$outcome <<- outcome
        rotate()
      },
      
      center.me = function () if (restrict.and.shift) {
        shift <<- mean(b.vector)
        b.vector <<- b.vector - shift
      },
       
      reinitialize = function (n.nodes=NULL,
        edge.list=NULL, node.names=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          sr.rows <<- row.list.maker(edge.list)
        }

        if (n.groups > n.nodes) {
          warning ("MMSBM: Resetting number of groups to one less than the number of nodes.")
          n.groups <<- n.nodes - 1
          b.vector <<- rep(0, n.groups*(n.groups+1)/2)
          
          membership.node <<- sample(n.groups, n.nodes, replace=TRUE)
          membership.edge <<- draw.MMSB.from.nodes(edge.list, membership.node)
          membership.alpha0 <<- membership.alpha0[1:n.groups]
        }
        
        
        if (ncol(membership.node) != n.nodes) {
          message ("Reinitializing MMSBM Membership Fractions")
          membership.node <<- rdirichlet.block (matrix(membership.alpha0, nrow=n.groups, ncol=n.nodes))
        }
        if (nrow(membership.edge) != nrow(edge.list)) {
          message ("Reinitializing MMSBM Edge Memberships")
          membership.edge <<- draw.MMSB.from.nodes(edge.list, membership.node)
        }
        rotate()
        if (!is.null(node.names)) {
          if (length(node.names) == .self$n.nodes) node.names <<- node.names
        } else node.names <<- as.character(1:.self$n.nodes)

      },

      pieces = function (include.name=FALSE) {
        out <- list (b.vector=b.vector,
                     membership.edge=membership.edge,
                     membership.node=membership.node)
        class(out) <- "MMSBMout"
        #if (include.name) out <- c("SBM", out)
        out
      },

      show = function () {
        message("b.vector:"); print(b.vector)
        message("membership.edge:"); print(t(membership.edge))
        message("membership.node:"); print(membership.node)
        #message("mult.factor:"); print(mult.factor)
      },
      plot = function (memb=membership.node, block=symBlock(b.vector), ...) {
        block.membership.plot (memb, block, node.labels=node.names, ...)
      },
      plot.network = function (color=outcome, ...) {
        netplot (edge.list, color, node.labels=node.names, ...)
      },
      


      value = function () {
        sbm.matrix <- symBlock(b.vector)
        #mult.factor*
          sbm.matrix[membership.edge[,1] +
                     dim(sbm.matrix)[1]*(membership.edge[,2]-1)]
      },
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {   #slightly slower.
        sbm.matrix <- symBlock(parameters[[1]])
        #parameters[[3]]*
        sbm.matrix[parameters[[2]][edges,1] +
                   dim(sbm.matrix)[1]*(parameters[[2]][edges,2]-1)]
      },


      
      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},

      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },


      random.start = function () {

        membership.node <<-
          rdirichlet.block (matrix(membership.alpha0, nrow=n.groups, ncol=n.nodes))
        
        membership.edge <<- draw.MMSB.from.nodes(edge.list, membership.node)
        
        b.vector <<- rnorm(n.groups*(n.groups+1)/2, 0, 0.5)
#        rotate()
        
      },
      

      rotate = function () {
        rotation <- MMSBM.ID.rotation(membership.node, n.groups)
        membership.edge <<- matrix(rotation[c(membership.edge)], ncol=2)
        membership.node <<- membership.node[rotation,]
        b.vector <<- SBM.rotate.bvector(b.vector, rotation)
      },
      
      draw = function (verbose=0, as.if.single=FALSE) {

        if (length(outcome) != nrow(edge.list)) stop ("MMSBM: outcome and edge.list have different lengths.")

        #Hold me!
        b.matrix <- symBlock(b.vector)
        b.block <- membership.edge
        b.node <- membership.node
        #if (verbose>1) print(b.memb)

        
        # draw edge and node memberships. We can do them simultaneously for each node.
        for (ii in sample(1:n.nodes)) {
          lefty <- sr.rows[[ii]][which(edge.list[sr.rows[[ii]],1] == ii)]
          righty <- sr.rows[[ii]][which(edge.list[sr.rows[[ii]],2] == ii)]
          log.pp.mat <- t(sapply(1:n.groups, function(gg) {
            b.block[lefty, 1] <- gg
            b.block[righty, 2] <- gg
            piece <- b.matrix[b.block[c(lefty,righty),1] +
                              (b.block[c(lefty,righty),2]-1)*n.groups]
            dnorm(outcome[c(lefty,righty)], piece, sqrt(residual.variance), log=TRUE)
            # + log(b.node[gg,ii])
          }))

          if (!as.if.single) {
            #print("nAIS")
            picks <- apply(log.pp.mat + log(b.node[,ii]), 2, function(cc) {
              cc <- cc - max(cc); sample(1:n.groups, 1, prob=exp(cc))
            })
            if (length(lefty)>0) b.block[lefty,1] <- picks[1:length(lefty)]
            if (length(righty)>0) b.block[righty,2] <- picks[length(lefty) + 1:length(righty)]
          
          #node membership
            counts <- sapply(1:n.groups, function(gg)
                             sum(membership.edge[,1]==gg & edge.list[,1]==ii) +
                             sum(membership.edge[,2]==gg & edge.list[,2]==ii))

            b.node[,ii] <- rdirichlet.one (counts + membership.alpha0)
            
          } else {
            #print("AIS")
            
            #assuming even prior odds on each group, but all memberships for a node are the same.
            cc <- apply(log.pp.mat, 1, sum); cc <- cc - max(cc)
            pick <- sample(1:n.groups, 1, prob=exp(cc))
            if (length(lefty)>0) b.block[lefty,1] <- pick
            if (length(righty)>0) b.block[righty,2] <- pick
            b.node[,ii] <- 1/n.groups
            
          }
            
        }
        
#        if (verbose>1) print(b.memb)
        membership.edge <<- b.block
        membership.node <<- b.node
       
        b.vector <<- sapply(1:length(b.vector), function(bb) {
          picks <- unique(c(which(membership.edge[,1]==group.pairs[bb,1] &
                                  membership.edge[,2]==group.pairs[bb,2]),
                            which(membership.edge[,1]==group.pairs[bb,2] &
                                  membership.edge[,2]==group.pairs[bb,1])))
          if (length(picks) > 0) {
            var.b <- 1/(length(picks)/residual.variance + 1/b.vector.v[bb])
            mean.b <- var.b*(sum(outcome[picks])/residual.variance + b.vector.m[bb]/b.vector.v[bb])
            output <- rnorm(1, mean.b, sqrt(var.b))
          } else output <- rnorm(1, 0, 0.5)
          output
        })


        if (restrict.and.shift) {center.me()}
#        rotate()
        
      },

      gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1,
        make.random.start=FALSE,
        as.if.single=FALSE) {
        out <- list()
        
        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw(as.if.single=as.if.single);
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
            if (report.interval > 0) if (index %% report.interval == 0) message("MMSBM ",index)
          } else if (round(index)==index) {
            if (report.interval > 0) if (index %% report.interval == 0) message("MMSBM burnin ",index)
          }
        }
        return(out)
      },
      
      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.summary = function (gibbs.out) {
        membs <- matrix(apply(sapply(gibbs.out, function(gg) gg$membership.node), 1, mean), ncol=n.nodes)
        colnames(membs) <- node.names
        
        bvec <- apply(sapply(gibbs.out, function(gg) gg$b.vector), 1, mean)
        return(list(membership.node=membs,
                    b.vector=bvec,
                    block=symBlock(bvec)))
      },
      print.gibbs.summary = function (gibbs.out) {
        get.sum <- gibbs.summary(gibbs.out)
        message ("Block membership mixes:")
        print (get.sum$membership.node)

        message ("Block value matrix:")
        print (get.sum$block)
        
        return(invisible(get.sum))
      },

      gibbs.plot = function (gibbs.out, ...) {
        get.sum <- gibbs.summary(gibbs.out)
        plot (get.sum$membership.node,
              get.sum$block,
              main = "MMSBM Summary from Gibbs Sampler", ...)
      },

      gibbs.node.colors = function (gibbs.out, colors=(1:n.groups) + 1) {
        rep("#DDDDFF", n.nodes)
        #get.sum <- gibbs.summary(gibbs.out)
        #return(colors[get.sum$modal.membership])
      }


      )
    )

