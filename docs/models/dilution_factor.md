(dilution_factor)=

# Dilution factor

Let's suppose we want to measure the transit depth of a planet orbiting the primary star in a binary system.
If the spatial resolution of our instrument is not good enough to completely separate the two stars (i.e., the projected pixel size on the plane of the sky is comparable with the angular separation of the two stars), the light of the target and the light of the secondary will fall into the same aperture. As a consequence, the planetary transit will result shallower, as the baseline flux will be given to the sum of the fluxes of the two stars. This is valid for any star entering the photometric aperture, not just in the case of a binary star.


If the magnitude of the target and contaminant tars in the passband of the instrument are known, it is possible to compute the flux ratio added by the contaminants and include it in the model. This factor is here called *dilution factor*.

In `PyORBIT`, the dilution factor is modeled as the ratio between the total flux of the *contaminant* stars over the flux of the *target* star. This facto is sometimes called *contamination ratio*.
Given $N$ contaminant stars whose flux is fully captured by the aperture around the target star, the dilution factor is given by this equation:

```{math}
:label: dilution_factor_fluxes
d = \frac{\sum_{i=1}^{N} F_i}{F_{\mathrm{target}}}
```

With this definition, the dilution factor is equal to zero in the absence of contaminants, equal to one if the flux of the contaminants is equal to the flux of the target, etc.. 
Be sure to always use the same bandpasses for the contaminant and reference stars.


```{warning}
The are several definition of the dilution factor in the literature. Remember to always check the definition before taking a value from the paper. 
```

If you are observing with TESS, and you are *certain* that all the lights of the contaminants is falling within the aperture, the dilution factor can be computed by using the TESS magnitudes $T$ from the [TESS input catalog](https://tess.mit.edu/science/tess-input-catalogue/):
```{math}
:label: dilution_Factor_magnitudes
d = \frac{\sum_{i=1}^{N} 10^{(-0.4 * T_{i})}}{10^{(-0.4 * T_{\mathrm{target}})}}
```

Generally speaking,  you should compute the fraction of light of the stellar Point Spread Function that falls into a given photometric aperture ($k$ parameter) for each star (target and contaminants), and then derive the dilution factor (and its associated error) following equation {numref}`dilution_factor_fluxes`, e.g., as done in [Mantovan et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.516.4432M/abstract).


### Prior on the dilution factor

It is not possible to determine the value of the dilution factor from the dataset alone. The dilution factor is indeed one of the main reasons why false positives as eclipsing binaries or transiting brown dwarfs may be mistaken as planets. If you don't include a prior on the dilution factor, likely you will get a strong degeneracy on the scaled planetary radius or the normalization factor, among other parameters.


### Combined use of normalization and dilution factors

In `PyORBIT`, the dilution factor will be added to the planetary transit model, where the baseline flux of the start (i.e., the out-of-transit flux) is assumed to be unitary. Hence, the complete *model* will have a baseline equal to $1. + d_{\mathrm{factor}}$.
If you flattened your lightcurve with a filtering algorithm or divided it by a value to bring its average value around one, the higher baseline model due to the presence of the dilution factor will not match anymore the baseline flux of your observations. This difference must be absorbed by including a normalization factor, which will multiply the model in order to match the data. For example, if you have a dilution factor equal to $d=0.25$, the baseline of the model will be equal to $1.25$ (in units of target flux), the normalization factor will be equal to $1./1.25=0.80$, so that the combination of transit model plus dilution, all multiplied by the normalization facto, will match the value of the flattened lightcurve (ideally centered around 1.).
In conclusion, **it is strongly advised to always use a normalization factor in combination with the dilution factor** 

## Model definition and requirements

- model_name: ``dilution_factor`` or ``local_dilution_factor``
- required common objects: ``dilution_factor``

The boundaries of the dilution factor are automatically set to `[0.0, 1.0]`. While the dilution factor must be positive, it may be necessary to expand the upper boundaries for extreme cases when the integrated flux of the contaminants is similar or higher to the one of the target star.

## Keywords

This model does not require any keyword.


## Examples

These examples are based on the configuration files available in the folder [`lightcurve_matern32`](https://github.com/LucaMalavolta/PyORBIT_examples/tree/main/lightcurves_matern32) of the [PyORBIT examples](https://github.com/LucaMalavolta/PyORBIT_examples) repository.  For the sake of readability, we omit parts of the configuration files not affected by changes, and replace them with `...` .

in this example, we model the `Asiago` dataset with a dilution factor (specific to the passband of the Asiago observations) and a local normalization factor (different for every dataset).
Although not necessary, I have enlarged the upper boundary of the dilution factor as an example.

```yaml
inputs:
  ...
  LCdata_ASIAGO:
    file:  datasets/ASIAGO_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_asiago
      - normalization_factor
      - dilution_factor_ASIAGO
  ...
common:
  ...
  dilution_factor_ASIAGO:
    model: dilution_factor
    priors:
      d_factor: ['Gaussian', 0.2157, 0.0056]
    boundaries:
      d_factor: [0.00, 2.00]
...
models:
  normalization_factor:
    model: local_normalization_factor
    boundaries:
      n_factor: [0.8, 1.2]
  ...
```

```{tip}
``dilution_factor`` is a common model, as such it must be defined in the `common` section, otherwise priors and boundaries will be ignored.
```

It is possible to include a dilution factor for each dataset, rather than for each instrument, although I fail to see a plausible reason to do so.

```yaml
inputs:
  ...
  LCdata_CROW1:
    file:  datasets/CROW1_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_crow
      - celerite2_matern32_crow
      - normalization_factor
      - dilution_factor_CROW
  LCdata_CROW2:
    file:  datasets/CROW2_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_crow
      - celerite2_matern32_crow
      - normalization_factor
      - dilution_factor_CROW
  LCdata_CROW3:
    file:  datasets/CROW3_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_crow
      - celerite2_matern32_crow
      - normalization_factor
      - dilution_factor_CROW
  ...
common:
  ...
models:
  normalization_factor:
    model: local_normalization_factor
    boundaries:
      n_factor: [0.8, 1.2]
  dilution_factor_CROW:
    model: local_dilution_factor
    priors:
      d_factor: ['Gaussian', 0.2157, 0.0056]
  ...
```

