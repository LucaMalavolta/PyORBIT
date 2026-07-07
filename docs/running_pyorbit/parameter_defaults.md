(parameter_defaults)=

# Boundaries, spaces, and priors

Parameter boundaries are a prerequisite of MCMC/NS samplers, and together with priors and spaces, they must be specified for each parameter. The code automatically selects the values that work well in most cases, so you may need to do only minor tuning.  This page explains how to change the default behaviour of `PyORBIT`, and then it collects the default parameter properties used by `PyORBIT` when the configuration file does not override them with `boundaries`, `spaces`, `priors`, or `fixed`. Keep in mind that `PyORBIT` will activate only the parameters  relevant to the model you selected for your data 


## General description and examples

For *every parameter* involved in a fit, it is possible to specify 
- the `boundaries` of the parameter spaces;
- the `spaces` over which to explore the parameter space;
- the `priors` for a given parameter;
- if a parameter has a `fixed` value.

Sometimes it is possible to specify a prior for a derived parameter as well.




### Spaces

A parameter can be explored in *Linear* (or *Natural*) space when it spans a small range, or in *Logarithmic* space when it spans multiple orders of magnitude. A smaller base for the logarithm provides smaller jumps between orders of magnitude, and may be useful to better map the posterior.
Given $\theta$ the sampling parameter (the one explored by the sampler) and $\gamma$ the physical parameters (the one entering the model), the following options are available:

| Keyword | Description |  Transformation | Inverse | Notes |
| :--- | :--- | :--- | :--- | :-- |
| `Linear` | Linear or Natural space | $\theta = \gamma$ | $\gamma = \theta$  | Natural space, 1:1 match |
| `Log_Natural`  | Natural logarithm |  $\theta = \ln{\gamma}$  |  $\gamma = {\rm e}^\theta $ |  |
| `Log_Base2`   | Logarithm in base 2 |  $\theta = \log_{2}{\gamma}$ | $\gamma = 2^{\theta}$ | | 
| `Log_Base10`  | Logarithm in base 10 | $\theta = \log_{10}{\gamma}$ | $\gamma = 10^{\theta}$ | | 
| `Sine_Angle`  | The sine of the angle parameter   | $\theta = \sin {\gamma} $ | $\gamma = \arcsin {\gamma} $  |  Uniform prior in $\sin (\gamma)$ |
| `Cosine_Angle`| The cosine of the angle parameter | $\theta = \cos {\gamma} $ | $\gamma = \arccos {\theta} $ |   Uniform prior in $\cos (\gamma)$ |


```{code-block} yaml
:lineno-start: 33
:emphasize-lines: 6-8
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      spaces:
        P: Linear
        K: Linear
      boundaries:
        P: [0.2750, 0.2850]
        K: [0.001, 20.0]
        Tc: [57744.00, 57744.10]
      priors:
        P: ['Gaussian', 0.280324956, 0.000000067]
        Tc: ['Gaussian', 57744.071508, 0.000103]

```


### Boundaries

By specifying the `boundaries` of a parameter, the sampler will not be allowed to explore values outside this range. 
Boundaries can be specified using the keyword `boundaries` within the corresponding model.


```{code-block} yaml
:lineno-start: 33
:emphasize-lines: 9-12
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      spaces:
        P: Linear
        K: Linear
      boundaries:
        P: [0.2750, 0.2850]
        K: [0.001, 20.0]
        Tc: [57744.00, 57744.10]
      priors:
        P: ['Gaussian', 0.280324956, 0.000000067]
        Tc: ['Gaussian', 57744.071508, 0.000103]
```



### Priors

The prior expresses our knowledge of the parameter before analysing the data. The code will automatically retrieve the boundaries for a given parameter if they are needed for the prior calculation (e.g., as in the Uniform prior). Due to the different treatment of priors in Nested Sampling algorithms, some priors are available only in `Linear` space when used in NS analysis. Those cases are reported in the *Spaces (NS)* column

| Keyword | Description |  Parameters | Spaces (NS) |
| :--- | :--- | :--- | :--- | 
| `Uniform` | Uniform prior (default choice) | | All |
| `Gaussian` | Gaussian prior | 1) center 2) scale | | All |
| `PositiveHalfGaussian` | Right side of a Half-Normal distribution  | 1) center \ 2) scale | | All |
| `NegativeHalfGaussian` | Left side of a Half-Normal distribution | 1) center 2) scale | | All |
| `HalfGaussian` | Alias for `PositiveHalfGaussian`| 1) center 2) scale | | All |
| `TruncatedJeffreys` | Jeffreys prior, normalized between the boundaries | |  | 
| `Jeffreys` | Alias for `TruncatedJeffreys` prior | | |
| `TruncatedModifiedJeffreys` | Modified Jeffreys prior, normalized between the boundaries | 1) $a$ (linear for $\theta < a$ |  | 
| `ModifiedJeffreys` | Alias for `TruncatedModifiedJeffreys` | 1) $a$ | |
| `WhiteNoisePrior` | Alias for `TruncatedModifiedJeffreys` | 1) $a$ | |
| `TruncatedRayleigh` | Rayleigh prior, normalized between the boundaries | 1) scale $\sigma$ | | 
| `BetaDistribution` | Beta distribution | 1) shape parameter $\alpha$ 2) shape parameter $\beta$ | |
| `ComplementaryGaussian` | Bigaussian distribution symmetric around 90° |  1) center 2) scale | |  |
| `SymmetricGaussian` | Bigaussian distribution symmetric around 0 |  1) center 2) scale | |  |

Example:


```{code-block} yaml
:lineno-start: 33
:emphasize-lines: 13-15
common:
  planets:
    b:
      orbit: circular
      use_time_inferior_conjunction: True
      spaces:
        P: Linear
        K: Linear
      boundaries:
        P: [0.2750, 0.2850]
        K: [0.001, 20.0]
        Tc: [57744.00, 57744.10]
      priors:
        P: ['Gaussian', 0.280324956, 0.000000067]
        Tc: ['Gaussian', 57744.071508, 0.000103]
```


```{important}
Boundaries and priors in the YAML file must be written in the physical (or
natural) parameter space. This is also true when the parameter is sampled in a
logarithmic space. In that case, boundaries must be strictly positive.
```


### Where to specify parameter's properties 

TO  BE DONE




## Dataset systematics


In the tables below:

- `Uniform []` means a uniform prior over the listed boundaries, with no extra
  prior hyperparameters.
- `None` in the `Fixed` column means that the parameter has no predefined value when fixed by the user or the algorithm. 
- `dynamic` means that the numerical boundary is computed from the dataset or
  from another model at initialisation time.
- Pattern rows such as `poly_c{0..9}` represent all parameters in the indicated
  range.

Each dataset will include one offset and one jitter parameter for each active flag in
the input file. The actual parameter names are `offset_0`, `offset_1`, ... and
`jitter_0`, `jitter_1`, ...

| Parameter | Boundaries | Space | Prior | Fixed | Notes |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `offset_{n}` | dynamic | `Linear` | `Uniform []` | `[0.0, 0.0]` | See subsection below. |
| `jitter_{n}` | dynamic | `Linear` | `Uniform []` | `[0.0, 0.0]` | See subsection below. |
| `jitter_{n}` for transit-time or transit-duration datasets | dynamic | `Logarithmic` | `Uniform []` | `[0.0, 0.0]` | Same boundary logic as `jitter_{n}`, but sampled logarithmically. |

### Offset boundaries

The boundaries associated with each offset will be as large as possible to accommodate improbable or anomalous situations, such as eccentric planets or massive long-period companions with partial RV coverage. 
 
Given `x` as the array with the time of the observations, `y` as the array with the data, and `e` as the corresponding error, the boundaries for the offset are determined in the following way:

```python
x_range = np.max(x) - np.min(x)
y_range = np.max(y) - np.min(y)
y_trend = np.abs(y_range/x_range)
y_diff = np.abs(np.mean(x)-Tref) * y_trend + 1000.
offset_boundaries = [np.min(y) - 10.*y_diff, np.max(y) + 10.*y_diff]
```
where $T_{\rm ref}$ is the reference time specified in the parame 

Broad boundaries do not represent a problem when running an MCMC sampler, given an adequate starting point. When running `pyDE+emcee`, for example, you must allow a sufficiently large number of generations for the `pyDE` code to properly explore the parameter space; typically, 50000 generations is a good starting value.


### Jitter boundaries

The jitter boundaries are computed dynamically by taking the maximum error in the array as the reference and multiplying or dividing by 100 to set the upper or lower limits, respectively. 

```{code} python
jitter_boundaries = [min(e)/100, 100 max(e)]
```

The lower boundary is intentionally set to a value other than zero to avoid problems for datasets where the jitter parameter is explored in logarithmic space. Since jitter is added in quadrature to the formal errors, a value much smaller than the error has the same effect as no jitter at all. Hence, in most cases it is not necessary to change the default lower limit to zero. 

### Changing the boundaries

Offset and jitter boundaries are determined using the whole dataset, without distinction by the flag specified in the 4<sup>th</sup> and 5<sup>th</sup>  columns (see [Prepare a dataset file](prepare_dataset) ) 

Let's take the case when a dataset is using two flags for jitter and two flags for offset, e.g., when putting together data collected with HARPS and HARPS-N, but you want to use a single covariance matrix for the stellar activity. You can specify different boundaries for each distinct jitter and offset parameters:

```{code-block} yaml
:lineno-start: 1
inputs:
  RVdata:
    file: datasets/K2-141_RV_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
      - tinygp_quasiperiodic
    boundaries:
      jitter_0: [ 0.00,   10.0]
      jitter_1: [ 0.00,   20.0]
      offset_0: [ -3500.0, -3300.0]
      offset_1: [ -4000.0, -3000.0]
```

```text
----- dataset:  RVdata
jitter_0      id:  11  s:Linear      b:[      0.0000,      10.0000]   p:Uniform   []
jitter_1      id:  12  s:Linear      b:[      0.0000,      20.0000]   p:Uniform   []
offset_0      id:  13  s:Linear      b:[  -3500.0000,   -3300.0000]   p:Uniform   []
offset_1      id:  14  s:Linear      b:[  -4000.0000,   -3000.0000]   p:Uniform   []
```

The code will automatically assign the If you don't specify the boundaries for one of the parameters, the automa

```{code-block} yaml
:lineno-start: 8
:emphasize-lines: 2,4
    boundaries:
      #jitter_0: [ 0.00,   10.0]
      jitter_1: [ 0.00,   20.0]
      #offset_0: [ -3500.0, -3300.0]
      offset_1: [ -4000.0, -3000.0]
```

```text
----- dataset:  RVdata
jitter_0      id:  11  s:Linear      b:[      0.0120,    1088.3783]   p:Uniform   []
jitter_1      id:  12  s:Linear      b:[      0.0000,      20.0000]   p:Uniform   []
offset_0      id:  13  s:Linear      b:[ -13430.6237,    6632.5000]   p:Uniform   []
offset_1      id:  14  s:Linear      b:[  -4000.0000,   -3000.0000]   p:Uniform   []
```

You can apply the same boundaries to all the jitter/offset parameters of a dataset by removing the numerical subscript: 

```{code-block} yaml
:lineno-start: 8
:emphasize-lines: 2,3
    boundaries:
      jitter: [ 0.00,   10.0]
      offset: [ -3500.0, -3300.0]
```


```text
----- dataset:  RVdata
jitter_0      id:  11  s:Linear      b:[      0.0000,      10.0000]   p:Uniform   []
jitter_1      id:  12  s:Linear      b:[      0.0000,      10.0000]   p:Uniform   []
offset_0      id:  13  s:Linear      b:[  -3500.0000,   -3300.0000]   p:Uniform   []
offset_1      id:  14  s:Linear      b:[  -3500.0000,   -3300.0000]   p:Uniform   []

```
### Changing priors and spaces

The same tricks above with `boundaries`  also apply to `spaces` and `priors`, 
for example, specifying generic priors as below:

```{code-block} yaml
:lineno-start: 8
:emphasize-lines: 5,6
    boundaries:
      jitter: [ 0.00,   10.0]
      offset: [ -3500.0, -3300.0]
    priors:
      jitter: ['HalfGaussian', 0.00, 10.0]
      offset: ['Gaussian', -3450.0, 10.0]
```

will affect all the flag-specific parameters of your dataset:


```text
----- dataset:  RVdata
jitter_0      id:  11  s:Linear      b:[      0.0000,      10.0000]   p:HalfGaussian   [ 0. 10.]
jitter_1      id:  12  s:Linear      b:[      0.0000,      10.0000]   p:HalfGaussian   [ 0. 10.]
offset_0      id:  13  s:Linear      b:[  -3500.0000,   -3300.0000]   p:Gaussian   [-3450.    10.]
offset_1      id:  14  s:Linear      b:[  -3500.0000,   -3300.0000]   p:Gaussian   [-3450.    10.]
```


:::{caution}
:name: a-tip-reference
You cannot mix generic and flag-specific names! For example, this configuration:

```{code-block} yaml
    boundaries:
      jitter: [ 0.00,   10.0]
      offset: [ -3500.0, -3300.0]
      offset_1: [ -4000.0, -3000.0]
      jitter_1: [0.00,   20.0]
```
will produce unexpected results, including ignored priors and wrong boundary assignments.
:::

## Planets

Common object: `common: planets`.

| Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- |
| `P` | `[0.4, 100000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `K` | `[0.001, 2000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `Tc` | `[0.0, 1000.0]` | `Linear` | `Uniform []` | `None` |
| `Tc_Tref` | `[0.0, 1000.0]` | `Linear` | `Uniform []` | `None` |
| `mean_long` | `[0.0, 360.0]` | `Linear` | `Uniform []` | `0.0` |
| `e_coso` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `e_sino` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `sre_coso` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `sre_sino` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `e` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `omega` | `[0.0, 360.0]` | `Linear` | `Uniform []` | `90.0` |
| `M_Me` | `[0.05, 1000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `M_Ms` | `[1e-09, 0.5]` | `Log_Base2` | `Uniform []` | `None` |
| `Me_Ms` | `[0.1, 2000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `i` | `[0.0, 180.0]` | `Linear` | `Uniform []` | `90.0` |
| `Omega` | `[0.0, 360.0]` | `Linear` | `Uniform []` | `180.0` |
| `R_Rs` | `[1e-05, 0.5]` | `Linear` | `Uniform []` | `0.05` |
| `a_Rs` | `[1e-05, 500.0]` | `Linear` | `Uniform []` | `1.0` |
| `b` | `[0.0, 2.0]` | `Linear` | `Uniform []` | `0.0` |
| `lambda` | `[-180.0, 180.0]` | `Linear` | `Uniform []` | `0.0` |
| `phase_amp` | `[0.0, 0.5]` | `Linear` | `Uniform []` | `0.0` |
| `delta_occ` | `[0.0, 0.5]` | `Linear` | `Uniform []` | `0.0` |
| `phase_off` | `[-180.0, 180.0]` | `Linear` | `Uniform []` | `0.0` |
| `albedo` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `redist` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `insol` | `[0.0, 1000000000.0]` | `Log_Base2` | `Uniform []` | `1.0` |
| `apo_center` | dynamic | `Linear` | `Uniform []` | `None` |
| `apo_timescale` | `[10.0, 100000.0]` | `Log_Base10` | `Uniform []` | `None` |

## Star Parameters

Common object: `common: star: star_parameters`.

| Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- |
| `radius` | `[0.0, 2.0]` | `Linear` | `Uniform []` | `1.0` |
| `mass` | `[0.0, 2.0]` | `Linear` | `Uniform []` | `1.0` |
| `density` | `[0.0, 5.0]` | `Linear` | `Uniform []` | `1.0` |
| `i_star` | `[0.0, 180.0]` | `Linear` | `Uniform []` | `90.0` |
| `cosi_star` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `1.0` |
| `v_sini` | `[0.0, 200.0]` | `Linear` | `Uniform []` | `1.6` |
| `rotation_period` | `[1.0, 1000.0]` | `Linear` | `Uniform []` | `27.0` |
| `activity_decay` | `[10.0, 10000.0]` | `Linear` | `Uniform []` | `1000.0` |
| `temperature` | `[2000.0, 11000.0]` | `Linear` | `Uniform []` | `5777.0` |
| `line_contrast` | `[0.0, 100.0]` | `Linear` | `Uniform []` | `50.0` |
| `line_fwhm` | `[0.0, 12.0]` | `Linear` | `Uniform []` | `6.0` |
| `rv_center` | `[-300.0, 300.0]` | `Linear` | `Uniform []` | `0.0` |
| `veq_star` | `[0.0, 70.0]` | `Linear` | `Uniform []` | `1.6` |
| `alpha_rotation` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `0.6` |
| `convective_c1` | `[0.0, 5.0]` | `Linear` | `Uniform []` | `0.0` |
| `convective_c2` | `[-5.0, 0.0]` | `Linear` | `Uniform []` | `0.0` |
| `convective_c3` | `[-5.0, 5.0]` | `Linear` | `Uniform []` | `0.0` |

## Stellar Activity

Common object: `common: activity`.

| Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- |
| `Prot` | `[1.0, 1000.0]` | `Linear` | `Uniform []` | `None` |
| `Pdec` | `[1.0, 1000.0]` | `Linear` | `Uniform []` | `None` |
| `Pcyc` | `[50.0, 10000.0]` | `Linear` | `Uniform []` | `None` |
| `Oamp` | `[0.0001, 2.0]` | `Log_Base2` | `Uniform []` | `None` |
| `Hamp` | `[1e-08, 1000000.0]` | `Linear` | `Uniform []` | `None` |
| `Camp` | `[1e-08, 1000000.0]` | `Linear` | `Uniform []` | `None` |
| `sho_scale` | `[1.0, 1000.0]` | `Linear` | `Uniform []` | `None` |
| `sho_decay` | `[1e-08, 1000000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `sho_sigma` | `[1e-08, 1000000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `grn_period` | `[1.0, 1000.0]` | `Linear` | `Uniform []` | `None` |
| `grn_sigma` | `[1e-08, 1000000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `rot_sigma` | `[1e-08, 1000000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `rot_fmix` | `[0.001, 1.0]` | `Linear` | `Uniform []` | `None` |
| `rot_Q0` | `[1e-08, 1000000.0]` | `Log_Base10` | `Uniform []` | `None` |
| `rot_deltaQ` | `[1e-08, 1000000.0]` | `Log_Base10` | `Uniform []` | `None` |
| `Vc` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `Vr` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `Lc` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `Bc` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `Br` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `rot_amp` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `con_amp` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `cos_amp` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `cos_der` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `cyc_amp` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `cyc_der` | `[-500.0, 500.0]` | `Linear` | `Uniform []` | `None` |
| `matern32_sigma` | `[1e-06, 1000000.0]` | `Log_Base10` | `Uniform []` | `None` |
| `matern32_scale` | `[0.001, 1000.0]` | `Log_Base10` | `Uniform []` | `None` |
| `matern32_multigp_sigma` | `[-10000.0, 1000.0]` | `Linear` | `Uniform []` | `None` |
| `matern32_multigp_sigma_deriv` | `[-10000.0, 1000.0]` | `Linear` | `Uniform []` | `None` |
| `sin_P` | `[1.0, 1000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `sin_K` | `[0.001, 2000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `sin_f` | `[1.0, 1000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `grn_k{0..9}_period` | `[1e-08, 100.0]` | `Linear` | `Uniform []` | `0.0` |
| `grn_k{0..9}_sigma` | `[1e-08, 100.0]` | `Log_Base2` | `Uniform []` | `None` |
| `osc_k{0..9}_period` | `[1e-08, 100.0]` | `Linear` | `Uniform []` | `0.0` |
| `osc_k{0..9}_sigma` | `[1e-08, 100.0]` | `Log_Base2` | `Uniform []` | `None` |
| `osc_k{0..9}_Q0` | `[1.0, 1000.0]` | `Log_Base10` | `Uniform []` | `None` |

## Limb Darkening

Common objects: `ld_linear`, `ld_quadratic`, `ld_square-root`,
`ld_logarithmic`, `ld_exponential`, `ld_power2`, and `ld_nonlinear`.

| Common object | Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `ld_linear` | `ld_c1` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `None` |
| two-coefficient laws | `ld_c1` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `None` |
| two-coefficient laws | `ld_c2` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `None` |
| two-coefficient laws | `ld_q1` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `None` |
| two-coefficient laws | `ld_q2` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `None` |
| `ld_nonlinear` | `ld_c1` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `None` |
| `ld_nonlinear` | `ld_c2` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `None` |
| `ld_nonlinear` | `ld_c3` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `None` |
| `ld_nonlinear` | `ld_c4` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `None` |

The two-coefficient laws are `ld_quadratic`, `ld_square-root`,
`ld_logarithmic`, `ld_exponential`, and `ld_power2`.

## Normalization, Dilution, Offset, and Jitter

| Common object | Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `normalization_factor` | `n_factor` | `[1e-06, 1000000.0]` | `Log_Base2` | `Uniform []` | `0.0` |
| `dilution_factor` | `d_factor` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `common_offset` | `offset` | dynamic | `Linear` | `Uniform []` | `0.0` |
| `common_jitter` | `jitter` | dynamic | `Linear` | `Uniform []` | `0.0` |

`common_offset` and `common_jitter` inherit their numerical boundaries from the
datasets that use them. If more than one dataset shares the same common
parameter, the boundary is expanded to include all relevant datasets.

## Polynomial and Detrending Models

Common objects: `polynomial_trend`, `detrending`, and `lightcurve_detrending`.

| Common object | Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `polynomial_trend` | `x_zero` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `polynomial_trend` | `x_offset` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `polynomial_trend` | `poly_factor` | `[-100.0, 100.0]` | `Linear` | `Uniform []` | `0.0` |
| `polynomial_trend` | `poly_c{0..9}` | `[-10.0, 10.0]` | `Linear` | `Uniform []` | `0.0` |
| `detrending` | `det_linear` | `[-10.0, 10.0]` | `Linear` | `Uniform []` | `0.0` |
| `detrending` | `det_poly` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `detrending` | `det_c{0..9}` | `[-100.0, 100.0]` | `Linear` | `Uniform []` | `0.0` |
| `detrending` | `det_m32_sigma` | `[1e-06, 1000000.0]` | `Log_Base10` | `Uniform []` | `None` |
| `detrending` | `det_m32_rho` | `[1e-06, 1000000.0]` | `Log_Base10` | `Uniform []` | `None` |
| `detrending` | `x_zero` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `lightcurve_detrending` | `coeff_linear` | `[-10.0, 10.0]` | `Linear` | `Uniform []` | `0.0` |
| `lightcurve_detrending` | `coeff_poly` | `[-10.0, 10.0]` | `Linear` | `Uniform []` | `0.0` |
| `lightcurve_detrending` | `coeff_c0` | `[-100000.0, 100000.0]` | `Linear` | `Uniform []` | `0.0` |
| `lightcurve_detrending` | `x_zero` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `lightcurve_detrending` | `poly_c{0..9}` | `[-1000000.0, 1000000.0]` | `Linear` | `Uniform []` | `0.0` |

## Correlation Models

Common objects: `correlation` and `complex_correlation`.

| Model or common object | Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `correlation` | `x_zero` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `correlation` | `corr_c{0..10}` | `[-100000.0, 1000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `complex_correlation` | `x_zero` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `complex_correlation` | `corr_c{0..10}` | `[-100000.0, 1000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `local_correlated_jitter` | `x_zero` | `[-100000.0, 100000.0]` | `Linear` | `Uniform []` | `None` |
| `local_correlated_jitter` | `c{1..order}` | `[0.0, 1000000.0]` | `Linear` | `Uniform []` | `None` |

The `local_correlated_jitter` coefficients are model-local rather than common
object parameters. Their highest order is set by the model keyword `order`.

## Harmonics and Sinusoids

Common objects: `harmonics` and `sinusoid`.

| Common object | Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `harmonics` | `P` | `[0.5, 1000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `harmonics` | `T0` | `[-100.0, 100.0]` | `Linear` | `Uniform []` | `0.0` |
| `harmonics` | `phase` | `[0.0, 360.0]` | `Linear` | `Uniform []` | `0.0` |
| `harmonics` | `amp_S{1..5}` | `[1e-06, 1000000.0]` | `Log_Base2` | `Uniform []` | `0.0` |
| `harmonics` | `amp_C{1..5}` | `[1e-06, 1000000.0]` | `Log_Base2` | `Uniform []` | `0.0` |
| `harmonics` | `pha_S{n}` | `[0.0, 360.0/n]` | `Linear` | `Uniform []` | `0.0` |
| `harmonics` | `pha_C{n}` | `[0.0, 360.0/n]` | `Linear` | `Uniform []` | `0.0` |
| `sinusoid` | `sine_period` | `[0.4, 100000.0]` | `Log_Base2` | `Uniform []` | `None` |
| `sinusoid` | `sine_amp` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `None` |
| `sinusoid` | `sine_phase` | `[0.0, 360.0]` | `Linear` | `Uniform []` | `0.0` |
| `sinusoid` | `sine_offset` | `[0.0, 360.0]` | `Linear` | `Uniform []` | `0.0` |
| `sinusoid` | `x_zero` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `sinusoid` | `x_offset` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `sinusoid` | `poly_factor` | `[-1000000000.0, 1000000000.0]` | `Linear` | `Uniform []` | `0.0` |
| `sinusoid` | `poly_c{0..9}` | `[-1000000.0, 1000000.0]` | `Linear` | `Uniform []` | `0.0` |

## CHEOPS and CCF Parameters

Common objects: `cheops_modelling` and `ccf_parameters`.

| Common object | Parameter | Boundaries | Space | Prior | Fixed |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `cheops_modelling` | `scale_factor` | `[-5.0, 20.0]` | `Linear` | `Uniform []` | `1.0` |
| `cheops_modelling` | `ramp` | `[-100.0, 100.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdt` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `d2fdt2` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdbg` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdcontam` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdsmear` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdx` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdy` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `d2fdx2` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `d2fdxdy` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `d2fdy2` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdsinphi` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdcosphi` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdcos2phi` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdsin2phi` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdcos3phi` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `cheops_modelling` | `dfdsin3phi` | `[-1.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `ccf_parameters` | `contrast_m` | `[-10.0, 10.0]` | `Linear` | `Uniform []` | `0.0` |
| `ccf_parameters` | `contrast_q` | `[0.0, 1.0]` | `Linear` | `Uniform []` | `0.0` |
| `ccf_parameters` | `fwhm_m` | `[-10.0, 10.0]` | `Linear` | `Uniform []` | `0.0` |
| `ccf_parameters` | `fwhm_q` | `[0.0, 70.0]` | `Linear` | `Uniform []` | `0.0` |
| `ccf_parameters` | `rv_offset` | `[-10.0, 10.0]` | `Linear` | `Uniform []` | `0.0` |

## Model Families

Most entries in the `models` section do not define new boundaries directly.
They select which common-object parameters are sampled and whether those
parameters are shared among datasets or local to a dataset.

| Model family | Parameter defaults used |
| :--- | :--- |
| `radial_velocities`, `rv_planets`, `planetary_velocities` | `planets`; optionally `star_parameters` when planet mass parametrizations are used. |
| `transit_times`, `Tc_planets` | `planets`; TTV variants can create dataset-local `Tc` parameters whose default boundaries are the time span of each transit subset unless overridden. |
| `batman_transit`, `pytransit_transit`, `pytransit_dynamical`, and their TTV/subset variants | `planets`, `star_parameters`, and a limb-darkening common object. |
| `batman_transit_eclipse_phasecurve`, `spiderman_thermal` | `planets`, plus eclipse and phase-curve parameters from the `planets` table. |
| `gp_quasiperiodic`, `gp_quasiperiodic_alternative`, `gp_quasiperiodic_derivative`, `gp_quasiperiodic_cosine`, `gp_pyaneti_quasiperiodic`, `gp_framework_quasiperiodic` | `activity`; optionally `star_parameters` for shared `rotation_period` or `activity_decay`. |
| `tinygp_*`, `celerite2_*`, `spleaf_*` Gaussian-process models | `activity`; optionally `star_parameters` for shared stellar activity hyperparameters. |
| Multidimensional GP models | `activity`; coefficients such as `rot_amp`, `con_amp`, `cos_amp`, `cyc_amp`, and Matern-3/2 multidimensional coefficients are often dataset-local. |
| `polynomial_trend`, `shared_polynomial_trend`, `local_polynomial_trend`, `subset_polynomial_trend` | `polynomial_trend`. |
| `detrending`, `full_detrending`, `polynomial_detrending`, `exponential_detrending`, `detrending_matern32` | `detrending`. |
| `correlation`, `complex_correlation` | `correlation` or `complex_correlation`. |
| `local_correlated_jitter` | The model-local parameters listed in [Correlation Models](#correlation-models). |
| `common_offset`, `common_jitter` | `common_offset` and `common_jitter`, with dynamic boundaries from the datasets. |
| `dilution_factor`, `local_dilution_factor` | `dilution_factor`. |
| `normalization_factor`, `local_normalization_factor`, `subset_normalization_factor` | `normalization_factor`. |
| `harmonics` | `harmonics`. |
| `sinusoid`, `local_sinusoid`, `sinusoid_common_period`, `sinusoid_polynomial_modulation` | `sinusoid`. |
| `cheops_detrending`, `cheops_factormodel` | `cheops_modelling`. |
| `spectral_rotation` and subset variants | `star_parameters` and, when requested, `ccf_parameters`. |
| `rossitermclaughlin_*` models | `planets`, `star_parameters`, limb darkening, and in some implementations `ccf_parameters`. |
