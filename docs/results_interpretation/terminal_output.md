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


```{code} bash
All the chains are longer than 50\*ACF, but some are shorter than 100\*ACF
PyORBIT should keep running for at least     34713 more steps to reach 100\*ACF
Suggested value for burnin:  2064
```


```{admonition} Don't rely too much on the ACF
From my personal experience, when the analysis of light curve is involved the ACF analysis wil suggest you to keep running the sampler even if corvengence has been reached.
```