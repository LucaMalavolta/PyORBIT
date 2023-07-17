(terminal_output)=

# Interpreting terminal output

`PyORBIT_results` will produce several files depending on the flags provided at
runtime. The terminal output is always present

```{danger} 
This page is based on the ouput of `PyORBIT` version 9. Version 10 provides some more advanced outputs that are not yet described here
```

## Terminal output

First thing you will get is a list of a few characteristics of your system. This
section may change depending on your computer configuration and the employed sampler.

```text
PyORBIT v9.1.0

Python version in use:
3.9.12 (main, Apr  5 2022, 06:56:58)
[GCC 7.5.0]
 LaTeX disabled by default

emcee version:  3.1.2

 Reference Time Tref: 59200.0

 Dimensions = 17
 Nwalkers = 68

 Steps: 50000
```

If you used the `-all` flag, you will get several log messages regarding the status
of plot preparation and plot printing. Those are detailed in the corresponding section.

## Confidence intervals

Confidence intervals and/or associated errors are computed using the 15.865th
and the 84.135th percentiles of the distribution. The median value of the
distribution (50th percentile) is then subtracted from
these two values. For clarity, these values are always followed by the string
`(15-84 p)`.
The following values:

```text
jitter_0             25.641528         -1.200384         1.372232 (15-84 p)
```

read as $\mathrm{jitter}_0 = 25.6_{-1.2}^{+1.4} ms^{-1}$ (in case of radial
velocities).

```{note}
The code will not try to format the output according to the significant figures of a measurement. Unit measurements and priors will not be displayed, at least for now.
```

The associated error is not provided for the starting point of the MCMC analysis
or the MAP values.


## Basic model selection

The code provides the *log-probability* function, with its two components (the
*log-priors* and the *log-likelihood*) explicitly reported right below.

```text
 LN posterior:  -319.565121      -2.851297     2.207220 (15-84 p)

 Median log_priors     = -177.66943313873546
 Median log_likelihood = -136.37527830557443
```

The code computes the Bayesian Information Criterion (BIC), the Akaike
Information Criterion (AIC) and the AIC with a correction for small sample sizes
(AIC). These values are computed using the median value of the log-probability / log-likelihood.
Formally, these three criteria should be computed using the *log-likelihood*,
but I've seen several cases where the *log-probability* is used instead. The
appropriate choice is left to the user.

```text
 Median BIC  (using likelihood) = 355.9106421289588
 Median AIC  (using likelihood) = 298.75055661114885
 Median AICc (using likelihood) = 299.37171702070515

 Median BIC  (using posterior)  = 711.2495084064296
 Median AIC  (using posterior)  = 654.0894228886198
 Median AICc (using posterior)  = 654.7105832981761
```

The *maximum a posteriori* (MAP) set of parameters is computed by picking the
sampling corresponding to the maximum value of the log-probability distribution.
As in the previous steps, different model selection criteria using either the
log-likelihood or the log-probability are reported here, with the difference
that now the corresponding *MAP* value is used rather than the corresponding *median* value.

```text
 MAP log_priors     = -177.66943313873546
 MAP log_likelihood = -136.28166552722797

 MAP BIC  (using likelihood) = 355.7234165722658
 MAP AIC  (using likelihood) = 298.56333105445594
 MAP AICc (using likelihood) = 299.18449146401224

 MAP BIC  (using posterior)  = 711.0622828497367
 MAP AIC  (using posterior)  = 653.9021973319268
 MAP AICc (using posterior)  = 654.5233577414831
```

Finally, the code provides a suggestion regarding the use of either AIC or AICc
following the standard definition

```text
 AIC suggested over AICs because NDATA (   600 ) > 40 * NDIM (    13 )
```

```{warning}
No check over the applicability of BIS/AIC/AICc is performed! Be sure that the underlying assumptions are valid in your specific case.
```

## Autocorrelation analysis

If the posteriors have been computed using an MCMC sampler, as `emcee`, an autocorrelation analysis
is performed automatically using the methods implemented in the sampler. For
more information, please check the [Autocorrelation analysis and
convergence
page](https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr) from
the `emcee` documentation.

`PyORBIT` will first compute the *integrated autocorrelation time* for each
variable using the built-in method in `emcee`. This method automatically check
if the chains are long enough for a reliable autocorrelation analysis.
If everything is fine, you'll get this message:

```text
The chains are more than 50 times longer than the ACF, the estimate can be trusted
```

If the chains are too short, the `tolerance` parameter (i.e., the minimum number
of autocorrelation times needed to trust the estimate) will be lowered to 20,
with a warning and a(likely overestimated) suggestion for the appropriate chain length.

```text
***** WARNING ******
The integrated autocorrelation time cannot be reliably estimated likely the chains are too short, and ACF analysis is not fully reliable
emcee.autocorr.integrated_time tolerance lowered to 20
If you still get a warning, you should drop these results entirely
Chains too short to apply convergence criteria
They should be at least 2800*nthin = 280000
 ```

The maximum value of the integrated time is then used as step length (or step
size) in the following analysis


```text
Computing the autocorrelation time of the chains
Reference thinning used in the analysis: 100
Step length used in the analysis: 13*nthin = 1300
```

The convergence criteria are clearly stated out, for your benefit:

```text
Convergence criteria: less than 1% variation in ACF after 50 times the integrated ACF
At least 50*ACF after convergence, 100*ACF would be ideal
Negative values: not converged yet
```

Finally, you'll get a (badly formatted) table to identify misbehaving variables,
and a suggestion regarding the number of steps required to satisfy the criterion
reported above and the (minimum) value for the burn-in.

| sampler parameter       | ACF | ACF*nthin  | converged at  |  nteps/ACF |
| :--------------        | :-------- | :--------  | :--------  | :-------- |
| RVdata_offset_0      |   8.441  |    844.1  |      764     |    117.6 |
| RVdata_jitter_0      |  12.740  |   1274.0  |     2064     |     76.9 |
| BISdata_offset_0     |   8.183  |    818.3  |      764     |    121.3 |
| BISdata_jitter_0     |   9.140  |    914.0  |      764     |    108.6 |
| ... | ... | ...| ... | ... |

Legend:

- sampler parameter: the parameter explored by the sampler (it may be different
  from the parameter of the model)
- ACF: the integrated autocorrelation time, computed on the thinned chains.
- ACF*nthin: the integrated autocorrelation time, in units of sampler steps
- converged at: number of sampler steps at which the chain

```text
All the chains are longer than 50\*ACF, but some are shorter than 100\*ACF
PyORBIT should keep running for at least     34713 more steps to reach 100\*ACF
Suggested value for burnin:  2064
```

```{admonition} Don't rely too much on the ACF
From my personal experience, when the analysis of light curve is involved the ACF analysis wil suggest you to keep running the sampler even if convergence has been reached.
```

Two rows of asterisks delimit the end of this section

## Confidence intervals of the parameters

Confidence intervals of the posteriors are provided 34.135th
percentiles from the median on the left and right sides, in addition to the
median as well, as specified in the [previous section](#confidence-intervals).

This section is divided in three groups:

- *Statistics on the posterior of the sampler parameters*: the direct output of the
  sampler. If the parameters have been parametrized (including exploration in
  logarithmic  space), the values displayed here
  may not be physically meaningful.
- *Statistics on the model parameters obtained from the posteriors samples*:
  confidence intervals of the model parameters (i.e., before
  re parametrization) as they should appear in the physical model. If no
  reparametrization has been applied to a given parameter, the *model* posterior
  will be identical to the *sampler* posterior.
- *Statistics on the derived parameters obtained from the posteriors samples*:
  confidence intervals of those parameters that are not directly involved in the
  optimization procedure, but are derived from the model parameters through an
  analytical transformation.

The scheme will be repeated two times: once to print the confidence intervals of
sampler/model/derived parameters, and again to print the corresponding MAP values.


Let's take an example involving the determination of a planetary radius though
light curve fitting:

### Statistics on the posterior of the sampler variables

```text
====================================================================================================
     Statistics on the posterior of the sampler variables
====================================================================================================

----- common model:  b
P                0     -2.157173      -0.000002      0.000002 (15-84 p) ([-2.251539, -2.058894])
Tc               1  59144.616223      -0.000286      0.000285 (15-84 p) ([59144.600000, 59144.630000])
b                2      0.483160      -0.041884      0.034602 (15-84 p) ([ 0.000000,  2.000000])
R_Rs             3      0.020650      -0.000464      0.000448 (15-84 p) ([ 0.000010,  0.500000])
sre_coso         4      0.016803      -0.199241      0.192704 (15-84 p) ([-1.000000,  1.000000])
sre_sino         5      0.059049      -0.212933      0.198307 (15-84 p) ([-1.000000,  1.000000])
```

When reporting the sampler results, the name of the parameter is followed by a
number indicating its position in the array of parameters passed to the
optimization algorithm. These numbers are for internal use only (they my change
from run to run), but they can be useful to check the dimensionality of the
problem. At the end of the line, the boundaries of the parameters are reported:
these boundaries can be assigned by the users (as in the case of *P* and *Tc*)
or internally defined (as for *b* and *R_Rs*).
<!---
I should include a page in the documentation where the theta_parameters dictionary is detailed
-->
In this example you can see that the orbital period *P* is negative: this is not
an error, as the parameter is being explored in logarithmic space (specifically,
Log2) and not in Natural base.
Note also that eccentricity $e$ and argument of pericenter $\omega$ are
parametrized as $\sqrt{e} \sin \omega$ (*sre_sino*) and $\sqrt{e} \cos \omega$
(*sre_coso*) following
[Eastman et al. 2013](https://ui.adsabs.harvard.edu/abs/2013PASP..125...83E/abstract).
**I should make priors and spaces more explicit at
output**

### Statistics on the model parameters obtained from the posteriors samples

```text
====================================================================================================
     Statistics on the physical parameters obtained from the posteriors samples
====================================================================================================

----- common model:  b
P                 2.241951e-01     -3.188571e-07     2.983407e-07 (15-84 p)
Tc                59144.616223         -0.000286         0.000285 (15-84 p)
b                     0.483160         -0.041884         0.034602 (15-84 p)
R_Rs                  0.020650         -0.000464         0.000448 (15-84 p)
e                     0.060706         -0.042746         0.066179 (15-84 p)
omega               122.589193        -86.415897       136.651806 (15-84 p)
```

The main difference with respect to the previous table is that the columns of
the parameter index and the parameter boundaries are not listed here, as not all
the parameters may have been directly involved in the optimization procedure.
The orbital period *P* is now expressed in the correct unit (days) as it has
been converted from logarithmic to natural space.
, while the
other parameters *b*, *Tc*, and *R_Rs* are identical as they did not go through
any transformation.
The *argument of pericenter* (*omega*) and the eccentricity *e* are now
explicitly reported, as they are the actual parameters of the physical model.

In the case of a circular orbit, *sre_sino* and *sre_coso* would not be listed
in the sampler parameters output, while the eccentricity *e* and the argument of
periastron *\omega* would be listed as *fixed parameters* and without confidence
interval, as the model still require these terms although they are not involved
in the optimization procedure.

```text
...
e                 0.000000e+00
omega                90.000000
...
```

### Statistics on the derived parameters obtained from the posteriors samples

```text
====================================================================================================
     Statistics on the derived parameters obtained from the posteriors samples
====================================================================================================

----- common model:  b
a_Rs                  2.150072         -0.021830         0.021366 (15-84 p)
i                    77.018467         -1.040808         1.206541 (15-84 p)
mean_long           102.186949         -0.370244         0.392672 (15-84 p)
R_Rj                  0.125422         -0.002961         0.002917 (15-84 p)
R_Re                  1.405843         -0.033195         0.032702 (15-84 p)
T_41                  0.031653         -0.000485         0.000543 (15-84 p)
T_32                  0.029880         -0.000534         0.000604 (15-84 p)
a_AU_(rho,R)          0.006241         -0.000080         0.000080 (15-84 p)
```

This table looks a lot like the *physical parameters* table, the main difference
with the difference that these values are either obtained though internal
transformation of the physical parameters, as in the case of the scaled
semi-major axis *a_Rs* or the orbital inclination *i*, or using some external
parameters, as the planetary radius in Jupiter radii *R_Rj*  and Earth radii
*R_Re*, both requiring a prior on the stellar radius which is however not
involved in the optimization procedure.

## Differences between MCMC and nested sampling outputs

The output described above applies to MCMC samplers. Nested sampling algorithms
will have slightly different output:

- Autocorrelation function analysis will not be performed.
- An extra step reporting the boundary of the parameters will be printed.
- A short summary with the confidence interval for the Bayesian evidence may be printed.
- Some additional plots may be produced according to the methods implemented in
  the specific sampler.

The additional outputs follow closely the examples reported in the
documentation of each sampler, so they will not be detailed here.

## Note on reparametrization and transformation of parameters

Median values and confidence intervals for model and derived parameters are
computed directly on the transformed posterior, rather than on the reported
values for of the sampler parameters. In other words, each set of parameters
from the posterior distribution is converted into the model/derived parameters,
generating a new distribution for each new parameter, then the median and the
confidence intervals are computed on newly generated distribution.
For this reason, the median of a model parameter may not correspond to the
direct conversion of the median of the sampler parameters. Using the example
shown before:

$$\sqrt{e} \sin{\omega} = 0.017 \\
  \sqrt{e} \cos{\omega} = 0.059$$
$$e_{\mathrm{sampler}} = (\sqrt{e} \sin{\omega})^2 + (\sqrt{e} \cos{\omega} )^2
  = 0.004$$

while the reported median value for the eccentricity in the model parameters
section is (correctly) $e = 0.061$.
