CEC
===

The R Package CEC performs clustering based on the cross–entropy clustering (CEC) method, which was recently developed with the use of information theory. The main advantage of CEC is that it combines the speed and simplicity of k-means with the ability to use various Gaussian mixture models and reduce unnecessary clusters.

In the basic use of this package the input dataset <tt>data</tt> and the initial number <tt>centers</tt> of clusters: <br />

```R
cec(x = ..., centers = ...)
```
have to be specified. Below, a simple session with <b>R</b> is presented, where the component
(waiting) of the Old Faithful dataset is split into two clusters.

```R
library("CEC")
attach(faithful)
cec <- cec(matrix(faithful$waiting), 2)
print(cec)
```

As its main outcome <b>CEC</b> returns data cluster membership <tt>cec\$cluster</tt>. The following parameters of 
clusters are obtained as well:
<ul>
<li> a list of means <tt>cec\$centers</tt>, </li>
<li> a list of covariances <tt>cec\$covariances.model</tt> </li>
<li> a list of probabilities <tt>cec\$probability</tt>. </li>
</ul>

Some additional information concerning the number of iterations, cost (energy) function and the number of clusters during the following iterations is also obtained.

Below, a session of \proglang{R} is presented which shows how to use the above parameters for plotting the data and the Gaussian models corresponding to the clusters.

```R
hist(faithful$waiting, prob = TRUE, main = "Histogram of Time between Old Faithful eruptions", xlab = "Minutes", ylim = c(0, 0.05));
for(i in c(1:2)){
  curve(cec$probability[i] * dnorm(x, mean = cec$centers[i], sd = sqrt(cec$covariances[[i]][1])), add = T, col = i + 1)  
}
```

As it was said, the discussed method, analogously to k-means, depends on the initial clusters memberships. Therefore, the initialization should be started a few times, which can be obtained with the use of parameter <tt>nstart</tt> 
```R
cec <- cec(x = ...,  centers = ..., nstart = ...)
```
The initial cluster membership function can be chosen by the use of <tt>centers.init</tt> either randomly, <tt>"random"</tt>, or with the method given by the k-means++ algorithm <tt>"kmeans++"</tt>. 

Two more parameters are important in the initialization. The first <tt>iter.max = 100</tt> equals the maximum number of iterations in one CEC start and the second  <tt>card.min = "5\%"</tt> is the percentage of the minimal size of each cluster. The second parameter specifies the minimal possible number points in each cluster (clusters which contains less points are removed). Since each cluster is described by a covariance matrix, the number of elements in the cluster must be larger than the dimension of the data.

One of the most important properties of the CEC algorithm is that it can be applied for various Gaussian models. Therefore, the <b>CEC</b> package includes the implementation of six Gaussian models, which can be specified by the parameter <tt>type</tt>.

General Gaussian distributions
===

The family containing all Gaussian distributions is considered first. 
The results of the general Gaussian CEC algorithm give similar results to those obtained  by the Gaussian Mixture Models. 
However, the authors' method does not use the EM (Expectation Maximization) approach for minimization but a simple iteration process (Hartigan method). Consequently,  larger datasets can be processed in shorter time.

The clustering will have the tendency to divide
the data into clusters in the shape of ellipses
(ellipsoids in higher dimensions). 
 
```R
library("CEC")
data("fourGaussians")
cec <- cec(fourGaussians, centers = 10, type = "all", nstart = 20)
plot(cec, xlim = c(0, 1), ylim = c(0, 1), asp = 1)
cec.plot.cost.function(cec)
```

Spherical CEC
===

The second family discussed contains spherical Gaussian distributions $\G_{(\cdot I)}$ which 
can be accessed by
```R
\code{cec(x = ..., centers = ..., type = "spherical")}
```
The original distribution will be estimated by spherical (radial) densities, which will result with splitting the data into circle-like clusters of arbitrary sizes (balls in higher dimensions). 

```R
cec <- cec(x = Tset, centers = 10, type = "spherical")
plot(cec, xlim = c(0, 1), ylim = c(0, 1), asp = 1)
```

Spherical CEC with fixed radius
===

The next model implemented in the <b>CEC</b> package is a spherical model with a fixed covariance: 
```R
cec(x = ..., centers = ..., type = "fixedr", param = ...)
```
Similarly to the general spherical model, the dataset will be divided into clusters resembling full circles, but with the radius determined by <tt>param</tt>.

```R
cec <- cec(x = Tset, centers = 10, type = "fixedr", param = 0.01)
plot(cec, xlim = c(0, 1), ylim = c(0, 1), asp = 1)
```

Diagonal CEC
===

The fourth model is based on diagonal Gaussian densities 
```R
cec(x = ..., centers = ..., type = "diagonal") 
```

In this case, the data will be described by ellipses for which the main semi-major axes are parallel to the axes of the coordinate system. 

```R
cec <- cec(x = Tset, centers = 10, type = "diagonal")
plot(cec, xlim = c(0, 1), ylim = c(0, 1), asp = 1)
```
 
Fixed covariance CEC
===

The next model contains Gaussians with an arbitrary fixed covariance matrix  e.g \newline
```R
cec(x = ..., centers = ..., type = "covariances", param = ...)
```
In this example <br />
<prev>
0.04  0
0     0.01
</prev>
is used, which means that the data is covered by fixed ellipses.

```R
cec <- cec(x = Tset, centers = 10, type = "covariance",  param = matrix(c(0.04, 0, 0, 0.01), 2))
plot(cec, xlim = c(0, 1), ylim = c(0, 1), asp = 1)
```

Fixed eigenvalue CEC
===

The last model is based on Gaussians with arbitrary fixed eigenvalues  
```R
cec(x = ..., centers = ..., type = "eigenvalues", param = ...) 
```
In this example <tt>lambda_1=0.01</tt>, <tt>lambda_2=0.001</tt> are used, which means that the set is covered by ellipses with fixed semi axes (which correspond to the eigenvalues). 

```R
cec <- cec(x = Tset, centers = 10, type = "eigenvalues", param=c(0.01, 0.001))
plot(cec, xlim = c(0, 1), ylim = c(0, 1), asp = 1)
```



