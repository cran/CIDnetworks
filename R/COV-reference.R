
#library(Rcpp); library(mvtnorm); library(msm); sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R");

# Covariates component: Reference Class

COVcid <-
  setRefClass(
    "COVcid",
    fields = list(
      covariates="matrix",
      coef.cov="numeric",
      coef.cov.m="numeric",
      coef.cov.V="matrix",
      coef.cov.P="matrix",

      cov.names="character",

      cov.block="matrix",

      node.names="character",
      n.nodes="numeric",
      outcome="numeric",
      edge.list="matrix",
      residual.variance="numeric",
      edge.list.rows="list"    #,
      ),

    methods=list(
      initialize = function (

        covariates=array(1, c(nrow(edge.list),1)),
        coef.cov=rep(0, ncol(covariates)),

        cov.names=as.character(colnames(covariates)),

        n.nodes=10,
        edge.list=make.edge.list(n.nodes),
        edge.list.rows=row.list.maker(edge.list),
        residual.variance=1,
        outcome=numeric(0),

        coef.cov.m=rep(0, ncol(covariates)),
        coef.cov.V=diag(1000000, rep(ncol(covariates), 2)),

        generate=FALSE

        ) {

        .self$n.nodes <<- n.nodes
        .self$edge.list <<- edge.list
        .self$edge.list.rows <<- edge.list.rows
        .self$residual.variance <<- residual.variance
        .self$node.names <<- as.character(1:.self$n.nodes)

        .self$covariates <<- as.matrix(covariates)
        .self$coef.cov <<- coef.cov
        if (ncol(.self$covariates) != length(coef.cov)) stop ("Covariate matrix columns: ", ncol(.self$covariates), ", coefficient vector length ", length(coef.cov))

        if (length(cov.names) != ncol(.self$covariates)) {
          warning ("Covariate matrix columns: ", ncol(.self$covariates), ", covariate names ", length(cov.names), ";replacing with defaults.")
          .self$cov.names <- paste0 ("cov", 1:ncol(.self$covariates))
        } else {
          .self$cov.names <<- cov.names
        }

        .self$coef.cov.m <<- coef.cov.m
        .self$coef.cov.V <<- as.matrix(coef.cov.V)

        .self$cov.block <<- t(covariates)%*%covariates

        .self$coef.cov.P <<- solve(coef.cov.V)

        if (generate) .self$generate() else .self$outcome <<- outcome

      },

      reinitialize = function (n.nodes=NULL, edge.list=NULL, node.names=NULL) {
        if (!is.null(n.nodes)) n.nodes <<- n.nodes
        if (!is.null(edge.list)) {
          edge.list <<- edge.list
          edge.list.rows <<- row.list.maker(edge.list)
        }

        if (!(class(covariates) %in% c("array", "matrix"))) stop ("COV: Covariates must be a matrix or array object.")

        if (nrow(covariates) == n.nodes & ncol(covariates) == n.nodes) {
          message("Detected a sociomatrix-style covariate array. Adjusting to match the edge list.")
          elements.1 <- edge.list[,1] + (edge.list[,2]-1)*n.nodes
          if (length(dim(covariates)) == 3)
            elements.1 <- c(outer(elements.1, (1:dim(covariates)[3] - 1)*n.nodes^2, "+"))
          covariates <<- matrix (covariates[elements.1], ncol=dim(covariates)[3])
        }

        #Now, we can correct!
        if (nrow(covariates) < nrow(edge.list)) {
          stop (paste0("The number of rows in the COV matrix, ",nrow(covariates)," is less than the number of edges, ",nrow(edge.list),"."))
        }
        if (nrow(covariates) > nrow(edge.list)) {
          stop (paste0("The number of rows in the COV matrix, ",nrow(covariates)," is greater than the number of edges, ",nrow(edge.list),"."))
        }
        if (!is.null(node.names)) {
          if (length(node.names) == .self$n.nodes) node.names <<- node.names
        } else node.names <<- as.character(1:.self$n.nodes)
      },

      pieces = function (include.name=FALSE) {
        out <- list (coef.cov=coef.cov)
        class(out) <- "COVout"
        out
      },

      show = function (show.cov=FALSE) {
        message("coef.cov:"); print(coef.cov)
        if (show.cov) {message("t(covariates):"); print(t(covariates))}
      },
      plot = function (coefs=coef.cov, names=cov.names, sd=NULL, interval=NULL, main = "Covariate Summary from Class Values", ...) {
        dotchart.coef (coefs, names, sd, interval, main=main, ...)
      },
      plot.network = function (color=outcome, ...) {
        image.netplot (edge.list, color, node.labels=node.names, ...)
      },


      value = function () {covariates%*%coef.cov},
      value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) {cbind(covariates[edges,])%*%parameters[[1]]},



      generate = function () {outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))},

      log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
        meanpart <- value.ext (parameters, edges)
        sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
      },



      random.start = function () {coef.cov <<- c(rmvnorm (1, coef.cov.m, coef.cov.V))},   #ncol(covariates)

      draw = function (verbose=0) {
        varblock <- solve(cov.block/residual.variance + coef.cov.P)
        meanblock <- varblock%*%(t(covariates)%*%outcome/residual.variance +
                                 coef.cov.P%*%coef.cov.m)
        coef.cov <<- c(rmvnorm(1, meanblock, varblock))
      },

      gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1, make.random.start=FALSE) {
        out <- list()
        if (make.random.start) random.start()
        for (kk in 1:(draws*thin+burnin)) {
          draw();
          index <- (kk-burnin)/thin
          if (kk > burnin & round(index)==index) {
            out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
            if (report.interval > 0) if (index %% report.interval == 0) message("COV ",index)
          }
        }
        return(out)
      },

      gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
        value.ext (gg)
      }),

      gibbs.summary = function (gibbs.out) {
        coef.cov.mat <- rbind(sapply(gibbs.out, function(gg) gg$coef.cov))
        ob1 <- data.frame(estimated.mean=round(apply(coef.cov.mat, 1, mean),3),
                          estimated.sd=round(apply(coef.cov.mat, 1, sd),3),
                          q2.5=round(apply(coef.cov.mat, 1, quantile, 0.025),3),
                          q97.5=round(apply(coef.cov.mat, 1, quantile, 0.975),3))
        rownames(ob1) <- cov.names
        return(ob1)
      },
      print.gibbs.summary = function (gibbs.sum) {
        message ("Coefficients:")
        print (gibbs.sum)
        return()
      },


      gibbs.plot = function (gibbs.out, ...) {
        get.sum <- gibbs.summary(gibbs.out)
        plot (get.sum[,1], interval = get.sum[,3:4], main = "Covariate Summary from Gibbs Sampler", ...)
      },

      gibbs.node.colors = function (gibbs.out) {
        rep("#DDDDFF", n.nodes)
      }

      )
    )


