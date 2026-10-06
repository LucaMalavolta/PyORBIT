(polynomial_normalization)=

# Polynomial normalization

A polynomial normalization describes a slowly varying *multiplicative*
light-curve baseline. It uses the same polynomial-trend implementations as
[Polynomial trends](../models/polynomial_trends.md), with
`normalization_model: True` in the model definition. For a transit light curve,
PyORBIT evaluates the model schematically as

```{math}
:label: polynomial_normalization_equation
F_{\mathrm{model}}(t) = P(t)\,[1 + F_{\mathrm{transit}}(t)],
\qquad
P(t) = \sum_{k=0}^{n} c_k
\left(\frac{t-x_0}{\Delta t}\right)^k .
```

Here `F_transit` is the deviation from the out-of-transit flux, `x_zero` is
the reference time $x_0$, and `time_interval` is $\Delta t$ in the same units
as the input times. The constant coefficient `poly_c0` sets the baseline at
`x_zero`; for a light curve normalized near one, keep it near one. The other
coefficients describe the baseline slope and curvature. Choose a low order
unless the data support a more complex baseline.

## Model definition and requirements

Use the existing polynomial-trend models and set `normalization_model: True`:

| Model name | Parameter scope | Photometric use |
| --- | --- | --- |
| `local_polynomial_trend` | One reference time and coefficient set per dataset. | Independent baselines for different visits, instruments or light curves. |
| `polynomial_trend` | One reference time and coefficient set shared by attached datasets. | One common baseline for datasets on the same time axis. |
| `subset_polynomial_trend` | One reference time and coefficient set per subset. | Separate baselines for visits stored in a single input file with subset flags. |

All variants use the `polynomial_trend` common object, created automatically
when it is absent from `common`. The `shared_polynomial_trend` variant also
supports normalization, but has extra amplitude and offset parameters; see
[Polynomial trends](../models/polynomial_trends.md) before using that variant.

The transit or eclipse model and the polynomial normalization must both be
listed under the photometric dataset's `models`. When a free `poly_c0` provides
the baseline, omit a separate free `normalization_factor` for the same dataset:
the two constant factors would be degenerate.

## Model parameters

| Name | Meaning |
| --- | --- |
| `poly_c0` | Multiplicative baseline at `x_zero`. |
| `poly_c1`, `poly_c2`, ... | Slope, curvature and higher-order terms. |
| `x_zero` | Reference time; PyORBIT fixes it automatically to `Tref` if it lies within the dataset, or otherwise to an average dataset time for a local polynomial. |
| `poly_subM_cN`, `x_zero_subM` | Coefficients and reference time for subset `M` with `subset_polynomial_trend`. |

Set boundaries for the coefficients in the model definition. Coefficients of
`local_polynomial_trend` are fitted independently for each dataset even when
they use the same model name. Dataset-specific boundaries can be placed in a
subsection named after the dataset. See [default polynomial parameter
properties](../running_pyorbit/parameter_defaults.md#polynomial-and-detrending-models)
for the broad built-in boundaries.

## Keywords

| Keyword | Photometric setting |
| --- | --- |
| `normalization_model` | Set to `True` to multiply the light-curve model by the polynomial. |
| `order` | Highest polynomial order; defaults to `1`. |
| `starting_order` | Set to `0` to include `poly_c0`. This is selected automatically by `normalization_model: True` unless overridden. |
| `time_interval` | Scales elapsed time before polynomial evaluation; defaults to `1.0`. For day-based times, `0.1` gives coefficients per 0.1 day, while `10.0` gives coefficients per 10 days. |
| `x_zero` | Optional explicit reference time in the same time system as the dataset. |

Do not set `exclude_zero_point: True` or `starting_order: 1` when the
polynomial must also determine the light-curve baseline: those settings remove
`poly_c0`. If the baseline is fixed by another model, they can be used to fit
only the varying part of the polynomial.

## Examples

### TESS transit light curve with a quadratic baseline

The [HD 189733 TESS quickstart
configuration](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/quickstart/HD189733_example03_TESSphotometry.yaml)
fits `HD189733_TESS_PyORBIT.dat` with `batman_transit`. Add the
`polynomial_normalization` entry below to that dataset and the `models`
section; retain the example's planet, star, limb-darkening and solver settings:

```yaml
inputs:
  LCdata_TESS:
    file: HD189733_TESS_PyORBIT.dat
    kind: Phot
    models:
      - lc_model
      - polynomial_normalization

models:
  lc_model:
    model: batman_transit
    limb_darkening: limb_darkening
    planets: [b]
  polynomial_normalization:
    model: local_polynomial_trend
    normalization_model: True
    order: 2
    time_interval: 365.25
    boundaries:
      poly_c0: [0.95, 1.05]
      poly_c1: [-0.05, 0.05]
      poly_c2: [-0.05, 0.05]
```

The TESS data in that example have day-based timestamps and flux near one.
The observations span roughly a year, so `time_interval: 365.25` makes the
slope and curvature coefficients refer to a one-year time scale. The
boundaries above are illustrative and should be
adjusted to the observed out-of-transit baseline.

### Independent baselines for ground-based light curves

The [joint light-curve and RV
example](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/lightcurves_matern32/LC_RV_model01.yaml)
uses a separate local normalization factor for its Asiago and CROW photometry.
To fit a linear baseline instead, replace `normalization_factor` in those
photometric datasets' model lists with `polynomial_normalization`, and add this
definition while retaining the example's transit, GP, planet and RV settings:

```yaml
inputs:
  LCdata_ASIAGO:
    file: datasets/ASIAGO_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_asiago
      - polynomial_normalization
  LCdata_CROW1:
    file: datasets/CROW1_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_crow
      - celerite2_matern32_crow
      - polynomial_normalization

models:
  polynomial_normalization:
    model: local_polynomial_trend
    normalization_model: True
    order: 1
    time_interval: 0.1
    boundaries:
      poly_c0: [0.8, 1.2]
      poly_c1: [-0.05, 0.05]
```

The source example also has TESS, CROW2 and CROW3 light curves. Apply the
same replacement in each one when using the shared model definition above.
Because the implementation is `local_polynomial_trend`, each dataset gets its
own `poly_c0`, `poly_c1` and `x_zero` despite sharing the definition. The
`time_interval` value expresses the slope per 0.1 day, appropriate for
night-long ground-based observations; tune the bounds to the actual data.
