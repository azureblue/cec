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

As its main outcome <b>CEC}</b> returns data cluster membership <tt>cec\$cluster</tt>. The following parameters of 
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

