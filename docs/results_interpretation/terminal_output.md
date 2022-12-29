(terminal_output)=

# Interpreting PyORBIT_Results output

PyORBIT_Results will produced several files depending on the flags provided at
runtime. The terminal output is always present

## Terminal output

First thing you will get is a list of a few characteristics of your system

```{code} bash

PyORBIT v9.0.20

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

### Associated errors {#associated_errors}

Associated errors are computed using the 15.865th and the 84.135th percentiles of
the distribution. The median value of the distribution (50th percentile) is then subtracted from
these two values. For clarity, these values are always followed by the string
`(15-84 p)`.
The follwing value:

```{code} bash
jitter_0             25.641528         -1.200384         1.372232 (15-84 p)
```

reads as $\mathrm{jitter}_0 = 25.6_{-1.2}^{+1.4} ms^{-1}$ in case of radial
velocities. 

```{note}
The code will not try to format the output according to the significant figures of a measurement. Unit measurements and priors will not be displayed, at least for now. 
```

The associated error is not provided for the starting point of the MCMC analysis
or the MAP values.


### Basic model selection 

The code provides the *log-probability* function, with its two components (the
*plog-riors* and the *log-likelihood*) explicitely reported right below.

```{code} bash
 LN posterior:  -319.565121      -2.851297     2.207220 (15-84 p)

 Median log_priors     = -177.66943313873546
 Median log_likelihood = -136.37527830557443
```

The code computes the Bayesian Information Criterion (BIC), the Akaike
Information Criterion (AIC) and the AIC with a correction for small sample sizes
(AIC). These values are computed using the medina value of the log-probability / log-likelihood.
Formally, these three criteria should be computed using the *log-likelihood*,
but I'vee seen several cases where the *log-probability* is used instead. The
appropriate choice is left to the user.

```{code} bash

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
that now the corresponding *MAP* value is used rather then the corresponding *median* value. 

```{code} bash

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

```{code} bash
 AIC suggested over AICs because NDATA (   600 ) > 40 * NDIM (    13 )
```

```{warning}
No check over the applicability of BIS/AIC/AICc is performed! Be sure that the underlying assumtpions are valid in your specific case. 
```

### Aurocorrelation analysis

If the posteriors have been computed using `emcee`, an autocorrelation analysis
is performed automatically using the methods implemented in the sampler. For
more information, please check the [Autocorrelation analysis and
convergence
page](https://emcee.readthedocs.io/en/stable/tutorials/autocorr/#autocorr) from
the `emcee` documentation.


`PyORBIT` will first compute the *integrated autocorrelation time* for each
variable using the built-in method in `emcee`. This method automatically check
if the chains are long enough for a reliable autocorrelation analysis.
If everything is fine, you'll get this message:

```{code} bash

The chains are more than 50 times longer than the ACF, the estimate can be trusted

```

If the chains are too short, the `tolerance` parameter (i.e., the minimum number
of autocorrelation times needed to trust the estimate) will be lowered to 20,
with a warning a warning and a(likely overestimated) suggestion for the appropriate chain length.

```{code} bash

***** WARNING ******
The integrated autocorrelation time cannot be reliably estimated likely the chains are too short, and ACF analysis is not fully reliable
emcee.autocorr.integrated_time tolerance lowered to 20
If you still get a warning, you should drop these results entirely
Chains too shoort to apply convergence criteria
They should be at least 2800*nthin = 280000
 ```

The maximum value of the integrated time is then used as step length (or step
size) in the following analysis


```{code} bash
Computing the autocorrelation time of the chains
Reference thinning used in the analysis: 100
Step length used in the analysis: 13*nthin = 1300
```

The convergence criteria are clearly stated out, for your benefit:

```{code} bash
Convergence criteria: less than 1% variation in ACF after 50 times the integrated ACF
At least 50*ACF after convergence, 100*ACF would be ideal
Negative values: not converged yet
```

Finally, you'll get a (badly formatted) table to identify misbehaving variables,
and a suggestion regarding the number of steps required to satisfy the criterion 
reported above and the (minimum) value for the burn-in.

| sample variable          | ACF | ACF*nthin  | converged at  |  nteps/ACF |
| :--------------        | :-------- | :--------  | :--------  | :-------- |
| RVdata_offset_0      |   8.441  |    844.1  |      764     |    117.6 |
| RVdata_jitter_0      |  12.740  |   1274.0  |     2064     |     76.9 |
| BISdata_offset_0     |   8.183  |    818.3  |      764     |    121.3 |
| BISdata_jitter_0     |   9.140  |    914.0  |      764     |    108.6 |
| ... | ... | ...| ... | ... |


- sample variable: the variable explored by the sampler (it may be different
  from the physical variables of the model)
- ACF: the integrated autocorrelation time, computed on the thinned chains. 
- ACF*nthin: the integrated autocorrelation time, in units of sampler steps
- converged at: number of sampler steps at which the chain 

```{code} bash
All the chains are longer than 50\*ACF, but some are shorter than 100\*ACF
PyORBIT should keep running for at least     34713 more steps to reach 100\*ACF
Suggested value for burnin:  2064
```

```{admonition} Don't rely too much on the ACF
From my personal experience, when the analysis of light curve is involved the ACF analysis wil suggest you to keep running the sampler even if corvengence has been reached.
```

Two rows of asteriks delimits the end of this section

### Confidence intervals of the parameters

 Confidence intervals of the posteriors are provided 34.135th
 percentiles from the median on the left and right sides, in addition to the
 median as well, as specified in the
 [previous section](#associated_errors)

This section is divided in three groups:
- Statistics on the posterior of the sampler variables: 
- Statistics on the physical parameters obtained from the posteriors samples 