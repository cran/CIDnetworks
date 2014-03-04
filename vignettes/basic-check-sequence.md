Checking the CIDnetworks Suite: Basic Modes
========================================================

A.C. Thomas, November 21

This vignette contains testing programs for the master CID class, with no subclasses attached -- the Gilbert model for random graphs, essentially. This is done to ensure that the underlying mechanics are correctly specified.



Basic Structure: Binary Data
----------------------------

The default mode for generating a CID model is binary data, not because it's the simplest case (Gaussian data is) but because it's the one most likely to appear in real network data.

The generative model is

$$
Z_{ij} \sim N(\mu, 1); Y_{ij} = \mathbb{I}(Z_{ij} > 0)
$$


```r
basic.b <- CID.generate(100, intercept = 1, generate = TRUE)  #class.outcome='binary' can also be used.
```

```
## Error: invalid assignment for reference class field 'edge.list', should be
## from class "matrix" or a subclass (was class "numeric")
```

```r
basic.b.pieces <- basic.b$pieces()
```

```
## Error: object 'basic.b' not found
```

```r
basic.b
```

```
## Error: object 'basic.b' not found
```


Let's see the edges:


```r
netplot(basic.b$edge.list, basic.b$outcome)
```

```
## Error: could not find function "netplot"
```


Now run the Gibbs sampler for this model directly.


```r
basic.b.gibbs <- basic.b$gibbs.full(draws = 200, report = 200, thin = 2, make.random.start = TRUE)
```

```
## Error: object 'basic.b' not found
```

```r
par(mfrow = c(2, 2))
basic.b$gibbs.plot(basic.b.gibbs)
```

```
## Error: object 'basic.b' not found
```



```r
basic.b.gibbs.unwind <- list.output.to.matrices(basic.b.gibbs)
```

```
## Error: could not find function "list.output.to.matrices"
```

```r
plot(basic.b.gibbs.unwind$intercept, main = "Intercept")
```

```
## Error: object 'basic.b.gibbs.unwind' not found
```

```r
abline(h = basic.b.pieces$intercept, col = 2, lwd = 3)
```

```
## Error: object 'basic.b.pieces' not found
```


The residual variance here is an artifact -- it always equals 1 in the binary case. Let's check it to be sure.


```r
plot(basic.b.gibbs.unwind$residual.variance, main = "Residual Variance")
```

```
## Error: object 'basic.b.gibbs.unwind' not found
```

```r
abline(h = basic.b.pieces$residual.variance, col = 2, lwd = 3)
```

```
## Error: object 'basic.b.pieces' not found
```


Finally, the log likelihood of the data at each iteration.


```r
plot(basic.b.gibbs.unwind$log.likelihood, main = "Log Likelihood")
```

```
## Error: object 'basic.b.gibbs.unwind' not found
```


Basic Structure: Gaussian data
------------------------------

Gaussian/Normal data: generate data with only an intercept and variance for 100 nodes. We must specify the class here. The data generating mechanism is

$$
Y_{ij} \sim N(\mu, \sigma^2)
$$


```r
basic <- CID.generate(100, intercept = 1, residual.variance = 2, generate = TRUE, 
    class = "gaussian")
```

```
## Error: invalid assignment for reference class field 'edge.list', should be
## from class "matrix" or a subclass (was class "numeric")
```

```r
basic
```

```
## Error: object 'basic' not found
```

```r
basic.pieces <- basic$pieces()
```

```
## Error: object 'basic' not found
```

```r
basic.gibbs <- basic$gibbs.full(draws = 200, report = 200, thin = 2, make.random.start = TRUE)
```

```
## Error: object 'basic' not found
```


Now, process the outcome from the Gibbs sampler.


```r
basic.gibbs.unwind <- list.output.to.matrices(basic.gibbs)
```

```
## Error: could not find function "list.output.to.matrices"
```

```r
basic.gibbs.value <- basic$gibbs.value(basic.gibbs)
```

```
## Error: object 'basic' not found
```

```r

plot(basic.gibbs.unwind$intercept, main = "Intercept")
```

```
## Error: object 'basic.gibbs.unwind' not found
```

```r
abline(h = basic.pieces$intercept, col = 2, lwd = 3)
```

```
## Error: object 'basic.pieces' not found
```

```r
plot(basic.gibbs.unwind$residual.variance, main = "Residual Variance")
```

```
## Error: object 'basic.gibbs.unwind' not found
```

```r
abline(h = basic.pieces$residual.variance, col = 2, lwd = 3)
```

```
## Error: object 'basic.pieces' not found
```

```r
plot(basic.gibbs.unwind$log.likelihood, main = "Log Likelihood")
```

```
## Error: object 'basic.gibbs.unwind' not found
```


Basic Structure: Ordinal Data
-----------------------------

Generate ordinal network edges with 4 categories using a Gaussian random variable and two additional positive cutoff values; the edge's category is determined by its position with respect to zero and the two cutoffs $c_1$ and $c_2$. In general, for $c \geq 3$ categories, the generative mechanism is

$$
Z_{ij} \sim N(\mu, 1); Y_{ij} = \mathbb{I}(Z_{ij} > 0) + \sum_{k=1}^{c-2} \mathbb{I}(Z_{ij} > c_k)
$$


```r
basic.o <- CID.generate(100, intercept = 1, generate = TRUE, class.outcome = "ordinal", 
    ordinal.count = 4)
```

```
## Error: invalid assignment for reference class field 'edge.list', should be
## from class "matrix" or a subclass (was class "numeric")
```

```r
basic.o
```

```
## Error: object 'basic.o' not found
```

```r
table(basic.o$outcome)
```

```
## Error: object 'basic.o' not found
```

```r
netplot(basic.o$edge.list, basic.o$outcome)
```

```
## Error: could not find function "netplot"
```

```r
basic.o.pieces <- basic.o$pieces()
```

```
## Error: object 'basic.o' not found
```

```r
basic.o.gibbs <- basic.o$gibbs.full(draws = 100, report = 100, burnin = 100, 
    thin = 5, make.random.start = TRUE)
```

```
## Error: object 'basic.o' not found
```



```r
basic.o.gibbs.unwind <- list.output.to.matrices(basic.o.gibbs)
```

```
## Error: could not find function "list.output.to.matrices"
```

```r
# par(mfrow=c(1,3))
plot(basic.o.gibbs.unwind$intercept, main = "Intercept")
```

```
## Error: object 'basic.o.gibbs.unwind' not found
```

```r
abline(h = basic.o.pieces$intercept, col = 2, lwd = 3)
```

```
## Error: object 'basic.o.pieces' not found
```

```r
plot(basic.o.gibbs.unwind$ordinal.cutoffs1, main = "Ordinal Cutoff 1|2")
```

```
## Error: object 'basic.o.gibbs.unwind' not found
```

```r
abline(h = basic.o.pieces$ordinal.cutoffs1, col = 2, lwd = 3)
```

```
## Error: object 'basic.o.pieces' not found
```

```r
plot(basic.o.gibbs.unwind$ordinal.cutoffs2, main = "Ordinal Cutoff 2|3")
```

```
## Error: object 'basic.o.gibbs.unwind' not found
```

```r
abline(h = basic.o.pieces$ordinal.cutoffs2, col = 2, lwd = 3)
```

```
## Error: object 'basic.o.pieces' not found
```

```r
plot(basic.o.gibbs.unwind$log.likelihood, main = "Log Likelihood")
```

```
## Error: object 'basic.o.gibbs.unwind' not found
```

