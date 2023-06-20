(correlated_datasets)=

# Correlated datasets

It is possible to incorporate a model for which the observed value $Y$ at epoch $i$ is a function of another (error-less) observable $Z$:


```{math}
:label: polynomial_correlation

Y_i = \sum_{p=0}^{d} c_p (Z_i-Z_\mathrm{ref})^p 

```

$Z_\mathrm{ref}$ is a reference value that can be either specified in the configuration file or automatically computed by `PyORBIT`. In the latter case, the average of all $Z$ is taken. 