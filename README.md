CEC
===

The R Package CEC performs clustering based on the cross–entropy clustering (CEC) method, which was recently developed with the use of information theory. The main advantage of CEC is that it combines the speed and simplicity of k-means with the ability to use various Gaussian mixture models and reduce unnecessary clusters.

In the basic use of this package the input dataset <tt>data</tt> and the initial number <tt>centers</tt> of clusters: <br />
```R
cec(x = ..., centers = ...)
```
have to be specified. Below, a simple session with <tt>R</tt> is presented, where the component
(waiting) of the Old Faithful dataset is split into two clusters.
```R
library("CEC")
attach(faithful)
cec <- cec(matrix(faithful$waiting), 2)
print(cec)
```
