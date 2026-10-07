(common_dilution_factor)=

# Dilution factor

The `dilution_factor` common object stores the contaminating flux used by a
photometric [dilution model](../photometry/dilution_factor.md). Its parameter
`d_factor` is the ratio of the flux from contaminating stars to the flux from
the target star, measured in the same passband:

```{math}
d_{\mathrm{factor}} = \frac{\sum_i F_{\mathrm{contaminant},i}}{F_{\mathrm{target}}}.
```

A value of zero means no contaminating light. Use a separate common object
for each passband with a different flux ratio; datasets in the same passband
can share one object.

```{note}
The dilution model adds `d_factor` to a transit model whose out-of-transit
target flux is one. When the observed light curve has been normalized to a
baseline near one, also fit a [normalization factor](../photometry/normalization_factor.md).
An external prior on dilution is useful because transit data alone cannot
usually distinguish it from a change in transit depth.
```

## Model definition and requirements

- container name: `common`
- common object name: `dilution_factor`, or a passband-specific name with `type: dilution_factor`
- source class: `CommonDilutionFactor`
- associated photometric model: `dilution_factor`

Declare the common object in `common` and include its name in the photometric
dataset's `models` list. PyORBIT then uses the same name for the photometric
model. For the dataset-specific alternative, see
[`local_dilution_factor`](../photometry/dilution_factor.md).

```yaml
common:
  dilution_factor_ASIAGO:
    type: dilution_factor
```

## Model parameters

| Name | Parameter | Unit |
| :--- | :-------- | :--- |
| `d_factor` | Contaminant-to-target flux ratio | unitless |

The default bounds are `[0.0, 1.0]`, with a uniform prior in linear space.
Increase the upper bound if the contaminating flux can exceed the target
flux. The source class declares `0.0` as the default fixed value, which is
used when the parameter is explicitly fixed to its default.

## Keywords

The common object does not require any model-specific keyword. When its name
differs from `dilution_factor`, select the source class as follows:

**type**
* accepted values: `dilution_factor` (required for a custom object name)
* selects the `CommonDilutionFactor` class. `model` and `kind` are accepted
  aliases for this class selector.

## Examples

The `lightcurves_matern32/LC_RV_model01.yaml` file in the
[PyORBIT examples](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/lightcurves_matern32/LC_RV_model01.yaml)
contains a commented dilution prior. The excerpt below activates that idea
for the Asiago passband. Other Asiago light curves can use the same common
object in their `models` list:

```yaml
inputs:
  LCdata_ASIAGO:
    file: datasets/ASIAGO_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_asiago
      - normalization_factor
      - dilution_factor_ASIAGO

common:
  dilution_factor_ASIAGO:
    type: dilution_factor
    boundaries:
      d_factor: [0.0, 2.0]
    priors:
      d_factor: ['Gaussian', 0.2157, 0.0056]

models:
  normalization_factor:
    model: local_normalization_factor
    boundaries:
      n_factor: [0.8, 1.2]
```

The upper bound of `2.0` is illustrative; choose bounds and a prior from
measurements in the observed passband. See the
[dilution model page](../photometry/dilution_factor.md) for the effect on the
light-curve baseline.
