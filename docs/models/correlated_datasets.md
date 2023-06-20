(correlated_datasets)=

# Correlated datasets

It is possible to incorporate a model for which the observed value $Y$ at epoch $i$ is a function of another (error-less) observable $Z$, which I will call *correlated variable*:


```{math}
:label: polynomial_correlation

Y_i = \sum_{p=0}^{d} c_p (Z_i-Z_\mathrm{ref})^p

```

$Z_\mathrm{ref}$ is a reference value that can be either specified in the configuration file or automatically computed by `PyORBIT`. In the latter case, the average of all $Z$ is taken.

## Model definition and requirements

- model_name: ``local_correlation`` or ``correlation``
- required common objects: ``correlation``

```{tip}
The common object ``correlation`` will be automatically grabbed by `PyORBIT`. Since there are no common parameters, there is no need to specify it in the ``common`` section
```

## Keywords

Model-wide keywords, with the default value in bold face.

**order**
* accepted values: any number up to 10 | **`1`**
* order (or degree) of the polynomial used to describe the correlation

**correlated variable**
* accepted values: any variable name  | **`corr`**
* specify which variable you want to use as correlated model. You need to declare it if you use a name different from `corr` for the label of your column in your dataset/ancillary file

**normalization_model**
* accepted values: `True` | **`False`**
* If `False`, the correlation model is *added* to the final model (which is the combination of all the other models). If `True`, the correlation model is *multiplied* to the final model (as in the case of a normalization factor).
* If set to `True`, they keyword `include_zero_point` will be set to True as well.

**use_median_xzero**
* accepted values: **`True`** | `False`
* if  `True `, automatically compute the median of the correlated variable, regardless of the value provided to the keyword `use_median_xzero`

**x_zero**
* accepted values: any value | **median of the correlated dataset**
* User value for $Z_\mathrm{ref}$ as defined in Equation {eq} `polynomial_correlation`. To be valid, the keyword `use_median_xzero` must be set to `False`.

**include_zero_point**
* accepted values: `True` | **`False`**
* The coefficient $c_0$ is automatically set to zero, in order to avoid a degeneracy between this parameter and the offset of the dataset. You can decide to fit for $c_0$ by setting this keyword to True, but be aware of the consequences
* Note: this keyword is set to `True` if `normalization_model` is set to True

**exclude_zero_point**
* accepted values: **`True`** | `False`
It does the opposite of the keyword `include_zero_point`. It must be set to `True` if you are using the `normalization_model` flag but you don't want to fit for $c_0$

**baseline_value**
* I'll fix this

```{warning}
Not all the keywords have been implemented in version <= 9.1.12, check out if the output is consistent with your expectations.
```
## Examples

To be completed

```yaml
inputs:
  RV_data:
    file: dataset_dat.dat
    ancillary: dataset_anc.dat
    kind: RV
    models:
      - radial_velocities
      - local_correlation
common:
  planets:
    b:
      orbit: circular
      boundaries:
        P: [2, 20]
        K: [0.01, 10.0]
  star:
    star_parameters:
      priors:
        radius: ['Gaussian', 1.00, 0.02]
        mass: ['Gaussian', 1.00, 0.02]
models:
  radial_velocities:
    planets:
      - b
  local_correlation:
    type: local_correlation
    order: 1
parameters:
  Tref: 0.00 #BJD-2450000; istante arbitrario, più o meno centrale
solver:
  pyde:
    ngen: 20000
    npop_mult: 4
  emcee:
    npop_mult: 4
    nsteps: 25000
    nburn: 10000
    thin: 100
  nested_sampling:
    nlive: 4000
  recenter_bounds: True
```
