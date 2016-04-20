CEC
===

The R Package CEC performs clustering based on the cross–entropy clustering (CEC) method, which was recently developed with the use of information theory. The main advantage of CEC is that it combines the speed and simplicity of k-means with the ability to use various Gaussian mixture models and reduce unnecessary clusters.

The CEC package is a part of CRAN repository and it can be installed by the following command:

```R
install.packages("CEC")
library("CEC")
```

The basic usage comes down to the function `cec` with two required arguments: input data (`x`) and the initial number of centers (`centers`):

```R
cec(x = ..., centers = ...)
```
Below, a simple session with **R** is presented, where the component
(waiting) of the Old Faithful dataset is split into two clusters:

```R
library("CEC")
attach(faithful)
cec <- cec(matrix(faithful$waiting), 2)
print(cec)
```

As the main result, CEC returns data cluster membership `cec$cluster`. The following parameters of 
clusters can be obtained as well:

- means (`cec$centers`)
- covariances (`cec$covariances.model`)
- probabilities (`cec$probability`)

Additional information concerning the number of iterations, cost (energy) function and the number of clusters at each iteration are also available.

Below, a session of **R** is presented which shows how to use the above parameters for plotting the data and the Gaussian models corresponding to the clusters.

```R
hist(faithful$waiting, prob = T, main = "Time between Old Faithful eruptions", xlab = "Minutes", 
  col = "lightgray", border = 0, ylim = c(0, 0.05))
  
for(i in c(1:2))
  curve(cec$probability[i] * dnorm(x, mean = cec$centers[i], sd = sqrt(cec$covariances[[i]][1])),
    add = T, col = i + 1, lwd = 2)  
```
![](https://azureblue.github.io/cec/static/old.png)

The CEC method, analogously to k-means, depends on the initial clusters memberships. Therefore, the initialization should be started a few times, which can be achieved using the `nstart` parameter.
```R
cec <- cec(x = ...,  centers = ..., nstart = ...)
```
The initial locations of the centers can be chosen either **randomly** or using the **k-means++** method and it's driven by the `centers.init` parameter which can take one of the two values: `"random"` or `"kmeans++"`.

In the context multiple starts, it is worth to mention the `iter.max` parameter which limits the number of iterations in each start.

An essential parameter, in the context of CEC method, is `card.min` which expresses minimal cluster size - the number of points, below which, the cluster is removed. Since each cluster is described by a covariance matrix, the number of elements in the cluster must be larger than the dimension of the data.

One of the most important properties of the CEC algorithm is that it can be applied to various Gaussian models. The CEC package includes the implementation of six Gaussian models, which can be specified by the parameter `type`.

Implemented Gaussian models (families)
--------------------------------------

### General Gaussian distributions
**`type = "all"`**

The general Gaussian CEC algorithm gives similar results to those obtained by the Gaussian Mixture Models. However, the CEC does not use the EM (Expectation Maximization) approach for minimization but a simple iteration process (Hartigan method). Consequently, larger datasets can be processed in shorter time.

The clustering will have a tendency to divide the data into clusters in the shape of ellipses (ellipsoids in higher dimensions). 
 
```R
data("fourGaussians")

cec <- cec(fourGaussians, centers = 10, type = "all", nstart = 20)

plot(cec, asp = 1)
cec.plot.cost.function(cec)
```
![](https://azureblue.github.io/cec/static/all.png)![](https://azureblue.github.io/cec/static/allcost.png)

### Spherical CEC 
**`type = "spherical"`**

The original distribution will be estimated by spherical (radial) densities, which will result in splitting the data into circle-like clusters of arbitrary sizes (balls in higher dimensions). 

```R
data("Tset")

cec <- cec(x = Tset, centers = 10, type = "spherical", nstart = 5)

plot(cec, asp = 1)
```
![](https://azureblue.github.io/cec/static/spherical.png)

### Spherical CEC with fixed radius
**`type = "fixedr"`**

Similarly to the general spherical model, the dataset will be divided into clusters resembling full circles, but with the radius determined by the `param` argument.

```R
data("Tset")

cec <- cec(x = Tset, centers = 10, type = "fixedr", param = 0.01, nstart = 20)

plot(cec, asp = 1)
```
![](https://azureblue.github.io/cec/static/fixedr.png)

### Diagonal CEC
**`type = "diagonal"`**

In this case, the data will be described by ellipses for which the main semi-major axes are parallel to the axes of the coordinate system. 

```R
data("Tset")

cec <- cec(x = Tset, centers = 10, type = "diagonal", nstart = 5)

plot(cec, asp = 1)
```
![](https://azureblue.github.io/cec/static/diagonal.png) 

### Fixed covariance CEC
**`type = "covariance"`**

This model contains Gaussians with an arbitrary fixed covariance matrix.

```R
data("Tset")

cec <- cec(x = Tset, centers = 10, card.min = '10%', type = "covariance",  
  param = matrix(c(0.04, 0, 0, 0.01), 2))

plot(cec, asp = 1)
```
![](https://azureblue.github.io/cec/static/cov.png)

In the above example, the following covariance matrix has been used, which results in covering the data by fixed shape ellipses:
```
0.04  0.00
0.00  0.01      
```

### Fixed eigenvalues CEC
**`type = "eigenvalues"`**

The last model is based on Gaussians with arbitrary fixed eigenvalues.

```R
data("Tset")

cec <- cec(x = Tset, centers = 10, type = "eigenvalues", param = c(0.01, 0.001), nstart = 5)

plot(cec, asp = 1)
```
![](https://azureblue.github.io/cec/static/eigen.png)

In the above example, two eigenvalues: **λ₁=0.01** and **λ₂=0.001** are used, which results in covering the data with ellipses having fixed semi axes (corresponding to the eigenvalues). 

A mix of the Gaussian models
----------------------------

One of the most powerful properties of the CEC algorithm is the possibility of mixing models. More precisely, the mixed models can be specified by giving a list of cluster types (and a list of parameters if needed).

```R
cec(x = ..., centers = ..., type = c("all", "diagonal", ...), param = ...).
```

```R
data("mixShapes")

cec <- cec(mixShapes, 7, iter.max = 3, 
  type = c("fixedr", "fixedr", "eigen", "eigen",  "eigen", "eigen", "eigen"),  
  param = list(350, 350, c(9000, 8), c(9000, 8), c(9000, 8), c(9000, 8), c(9000, 8)), nstart = 500)

plot(cec, asp = 1)
```
![](https://azureblue.github.io/cec/static/mix.png)

