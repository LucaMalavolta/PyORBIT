(common_limb_darkening)=

# Limb darkening

The `limb_darkening` common object holds the coefficients of a stellar
intensity law. Transit and Rossiter-McLaughlin models can share these
coefficients when they use the same passband. Define separate common objects
for passbands that need different coefficients.

```{note}
Only the coefficients required by the chosen law need to be configured. The
selected model determines which common parameters are sampled; unspecified
bounds and priors come from the source class.
```

## Model definition and requirements

- container name: `star`
- common object name: `limb_darkening` by default, or a descriptive passband name
- source classes: the `LimbDarkening_*` classes in `pyorbit/common/limb_darkening.py`

Each object is nested under `common.star`. Select a law with `model`, `type`,
or `kind`, then refer to the object by name in a transit model's
`limb_darkening` keyword. The available class selectors are:

| Selector | Law | Coefficients |
| :------- | :-- | :----------- |
| `ld_linear` | linear | one |
| `ld_quadratic` | quadratic | two |
| `ld_square-root` | square root | two |
| `ld_logarithmic` | logarithmic | two |
| `ld_exponential` | exponential | two |
| `ld_power2` | power two | two |
| `ld_nonlinear` | nonlinear | four |

```yaml
common:
  star:
    limb_darkening_TESS:
      model: ld_quadratic
models:
  lc_model_TESS:
    model: batman_transit
    planets: [b]
    limb_darkening: limb_darkening_TESS
```

## Model parameters

| Name | Parameter | Unit |
| :--- | :-------- | :--- |
| `ld_c1` | First coefficient of the chosen intensity law | unitless |
| `ld_c2` | Second coefficient, for two- and four-coefficient laws | unitless |
| `ld_c3` | Third coefficient of the nonlinear law | unitless |
| `ld_c4` | Fourth coefficient of the nonlinear law | unitless |
| `ld_q1` | First transformed coefficient in the Kipping parametrization | unitless |
| `ld_q2` | Second transformed coefficient in the Kipping parametrization | unitless |

The default bounds are `[0.0, 1.0]` for all listed parameters except
`ld_c2` in a two-coefficient law, whose bounds are `[-1.0, 1.0]`.
For a quadratic law with `parametrization: Kipping`, PyORBIT samples `ld_q1`
and `ld_q2` and derives the physical coefficients as
`ld_c1 = 2 sqrt(ld_q1) ld_q2` and
`ld_c2 = sqrt(ld_q1) (1 - 2 ld_q2)`.

## Keywords

The default keyword value is highlighted in boldface.

**model**
* accepted values: `ld_linear` | `ld_quadratic` | `ld_square-root` |
  `ld_logarithmic` | `ld_exponential` | `ld_power2` | `ld_nonlinear`
* selects the limb-darkening law. This selector is required when the common
  object has a descriptive name such as `limb_darkening_TESS`; `type` and
  `kind` are accepted aliases.

**parametrization**
* accepted values: **`Standard`** | `Kipping`
* uses the coefficients `ld_c1` and `ld_c2` directly with `Standard`, or
  samples `ld_q1` and `ld_q2` for a two-coefficient law with `Kipping`.
  The Kipping transformation above is commonly used with the quadratic law.

## Examples

The [single-band TESS example](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/quickstart/HD189733_example03_TESSphotometry.yaml)
uses a quadratic law with Kipping coefficients:

```yaml
common:
  star:
    limb_darkening:
      model: ld_quadratic
      parametrization: Kipping
models:
  lc_model:
    model: batman_transit
    planets: [b]
    limb_darkening: limb_darkening
```

The [multiband example](https://github.com/LucaMalavolta/PyORBIT_examples/blob/main/quickstart/HD189733_example04_multiband_photometry.yaml)
uses separate TESS and Cousins $I_c$ objects. It supplies priors on the
standard coefficients for the latter:

```yaml
common:
  star:
    limb_darkening_TESS:
      model: ld_quadratic
      parametrization: Kipping
    limb_darkening_Ic_Cousins:
      type: ld_quadratic
      priors:
        ld_c1: ['Gaussian', 0.45, 0.05]
        ld_c2: ['Gaussian', 0.13, 0.05]
models:
  lc_model_TESS:
    model: pytransit_transit
    planets: [b]
    limb_darkening: limb_darkening_TESS
  lc_model_Ic_Cousins:
    model: batman_transit
    planets: [b]
    limb_darkening: limb_darkening_Ic_Cousins
```
