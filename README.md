
![PyORBIT logo](docs/_static/PyORBIT_logo_transp.png?raw=true)

### A code for exoplanet orbital parameters and stellar activity.
***
## `PyORBIT` version 10 by Luca Malavolta - July 2024


## Documentation update

I have extensively updated the documentation and  several examples are now available in a new dedicated repository, check them out:
- [`PyORBIT` 10 Documentation](https://pyorbit.readthedocs.io/)
- [`PyORBIT` examples repository](https://github.com/LucaMalavolta/PyORBIT_examples)

For nostalgic people, `PyORBIT` 8 and 9 are available as branches of the main repository: [legacy_version8](https://github.com/LucaMalavolta/PyORBIT/tree/legacy_version8) and [legacy_version9](https://github.com/LucaMalavolta/PyORBIT/tree/legacy_version8) respectively.

## Updates on version 10

```{warning}
Starting from version 10.3, `PyORBIT` has been upgraded to support `tinygp` (version 0.3.0), which in turns requires Python **3.10** to work properly.
If you are using `PyORBIT` \=> 10.3, follow the installation instructions to create an environment with Python 3.10
```

- *Improved speed*
After several failed attempts, I finally managed to apply the advice from the [emcee parallelization page](https://emcee.readthedocs.io/en/stable/tutorials/parallel/) to the rather complex structure of PyORBIT. The speed-up is noticeable for large datasets (e.g., photometry).

- *Rossiter McLaughlin*
Rossiter McLaughlin effect can now be precisely modelled using the CCF simulation approach employed in [Covino et al. 2013](https://ui.adsabs.harvard.edu/abs/2013A%26A...554A..28C/abstract).
When the rotation period is known - together with the stellar radius - the stellar inclination can be derived avoiding the bias reported by [Masuda & Winn 2020](https://ui.adsabs.harvard.edu/abs/2020AJ....159...81M/abstract).
This model has been successfully employed in [Mantovan et al. 2024b](https://ui.adsabs.harvard.edu/abs/2024A%26A...684L..17M/abstract) for the characterization of TOI-5398b.

- *Multidimensional Gaussian Process*
The model was introduced a few years back, but now it can finally take advantage of improved parallelization.
Recent examples of multidimensional Gaussian Processes through `PyORBIT` can be found in [Nardiello et al. 2022](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract) and [Mantovan et al. 2024a](https://ui.adsabs.harvard.edu/abs/2024A%26A...682A.129M/abstract).

- *tinyGP for better performances*
Working on both *classic* and _multidimensional_ Gaussian Processes, although the former is showing some JAX problems when producing the output results.

```{text}
**No back-compatibility**
Version 10 is not compatible with the results obtained with version 9.
If you have been using the development version of V10, you may run into incompatibility issues as well.
```

## Updates on version 9

* Added celerite2 SHO term

* Improved output: `spaces`, `bounds`, and `priors` explicitely written out at runtime

* Minor changes:
  * `dataset_variables` and `common_variables` changed to `dataset_parameters`
    and `common_parameters`
  * New models added
  * All instances named `variables` renamed to `parameters` for internal
    consistency. Backward compatibility ensured by new method in `model_container_abstract`

* `PyORBIT` will use the correct parametrization when a prior is assigned to a parameter and nested sampling is used (NS do not allow to assign priors to derived parameters).

* Many new models are now available
  * Multidimensional Gaussian Process (also known as GP Framework, [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract), [Barragán et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..866B/abstract)).
  * New kernels for Gaussian Process regression: Quasi-Periodic with Cosine [(Perger et al. 2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..58P/abstract), Quasi-Periodic with first derivative.
  * CHEOPS detrending (similar to FactorModel as in [pycheops](https://github.com/pmaxted/pycheops)).
  * Phase curve and secondary eclipse modelling with `batman` [(Kreidberg 2015)](https://ui.adsabs.harvard.edu/abs/2015PASP..127.1161K/) or `spiderman` [(Louden et al. 2016)](https://ui.adsabs.harvard.edu/abs/2018MNRAS.477.2613L/abstract).
  * Rossiter-McLaughlin through analytical formulation by [Ohta et al. 2005](https://ui.adsabs.harvard.edu/abs/2005ApJ...622.1118O/abstract).
  * Celerite/Celerite2 _standard_ models, now better organized.
  * Fit of individual time of transits (with automatic definition of boundaries) for TTV analysis.
  * experimenting with `tinygp`

* More options for parameters space exploration: `Linear`, `Log_Natural`, `Log_Base2`, `Log_Base10`.

* Some variables have been renamed, to improve clarity of results:

  | Definition  | PyORBIT 8.x | PyORBIT 9.x |
  | ----------- | ----------- | ----------- |
  | Mean Longitude | f | mean_long |
  | Scaled planetary radius | R  | R_Rs |
  | Scaled semi-major axis | a | a_Rs |
  | Planetary mass in Earth masses | M | M_Me |
  | Stellar density | rho | density |

  Also, `batman_ld_quadratic` and `pytransit_ld_quadratic` have been merged into `ld_quadratic`.

  Note: the Mean Longitude is defined assuming the longitude of the ascending node $\Omega$ equal to zero, thus corresponding to the angle defined in section 4.3 of  [(Ford 2006)](https://ui.adsabs.harvard.edu/abs/2006ApJ...642..505F/abstract) and simply called _phase_ in [Malavolta et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016A%26A...588A.118M/abstract)

* Overall reorganization of the code

**Warning** Loss of backward-compatibility

  You cannot analyzes results obtain with previous versions (< 9) of PyORBIT. No worries, the old version is still available in the ```legacy``` branch, you can download it from the Github page or switch to it through the terminal:

  ```bash
  git checkout legacy
  ```
  To switch back to the current version, just execute:

  ```bash
  git checkout main
  ```

**Working on it**

  * Rossiter-McLaughlin through `starry` [(Luger et al. 2019)](https://ui.adsabs.harvard.edu/abs/2019AJ....157...64L/abstract)
  * Multi-component GP for light curves (Rotation + Granulation, Rotation + Granulation + Oscillations), following [Barros et al. 2020](https://ui.adsabs.harvard.edu/abs/2020A%26A...634A..75B/abstract).

**Samplers**

Bayesian evidence estimation can now be performed with:

  * `dynesty` by Josh Speagle [dynesty documentation](https://dynesty.readthedocs.io/en/latest/), [code repository](https://github.com/joshspeagle/dynesty/)
  * `UltraNest` by Johannes Buchner [UltraNest documentation](https://johannesbuchner.github.io/UltraNest/), [code repository](https://github.com/JohannesBuchner/UltraNest/)
  * Multinest and Polychord support is still there, but not supported anymore

  Just substitute "emcee" with "dynesty" or "ultranest" to when running the code.
  *Warning* MCMC and Nested Sampling handle prior in a radically different way, as such it is not possible to directly translate some priors from one sampler to another

**Documentation**
  Some incomplete documentation is available [here](http://pyorbit.readthedocs.io/). For any doubt, feel free to open an issue on GitHub - emails tend to be rather ineffective lately - I'll be happy to work out together any problem that may arise during installation or usage of this software.

  `PyORBIT` handles several kinds of datasets, such as radial velocity (RV), activity indexes, and photometry, to simultaneously characterize the orbital parameters of exoplanets and the noise induced by the activity of the host star. RV computation is performed using either non-interacting Kepler orbits or n-body integration. Stellar activity can be modeled either with sinusoids at the rotational period and its harmonics or Gaussian process. Offsets and systematics in measurements from several instruments can be modeled as well. Thanks to the modular approach, new methods for stellar activity modeling or parameter estimation can be easily incorporated into the code.

**Models**
Any of these models can be applied to a dataset. The user can choose which models should be used for each dataset.
- `Gaussian Processes` for RV or photometry (shared or independent hyperparameter)
- Transits, eclipses, phase curves
- `Polynomial trends` with user-defined order
- `Correlation` with activity indexes (or any other dataset)
- `Sinusoids` (independent or shared amplitudes and periods)
- `Harmonics` to test old results in the literature

**Priors**
These priors can be applied to any of the parameters (it's up to the user to choice the appropriate ones):

- `Uniform` default prior for all the parameters
- `Gaussian`
- `Jeffreys`
- `Modified Jeffreys`
- `Truncated Rayleigh`
- `WhiteNoisePrior`
- `BetaDistribution`

`Jeffreys` and `Modified Jeffreys` priors are actually `Truncated Jeffreys` and `Truncated Modified Jeffreys`, with truncation defined by the boundaries of the parameter space.

**Parameter exploration**
The user can choose between `Linear` and `Logarithmic`. Note that in the second case the parameter space is transformed into base-2 logarithm.

Most of the information can be found in [Malavolta et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016A%26A...588A.118M/abstract) and [Malavolta et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018AJ....155..107M/abstract).

**Papers using PyORBIT**

- Orbital obliquity of the young planet TOI-5398 b and the evolutionary history of the system [Mantovan et al. 2024b](https://ui.adsabs.harvard.edu/abs/2024A%26A...684L..17M/abstract)
- Gliese 12 b, a temperate Earth-sized planet at 12 parsecs discovered with TESS and CHEOPS [Dholakia & Palethorpe et al. 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.531.1276D/abstract)
- The GAPS Programme at TNG. LIII. New insights on the peculiar XO-2 system [Ruggieri et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...684A.116R/abstract)
- The GAPS programme at TNG. L. TOI-4515 b: An eccentric warm Jupiter orbiting a 1.2 Gyr-old G-star [Carleo et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...682A.135C/abstract)
- GAPS XLIX. TOI-5398, the youngest compact multi-planet system composed of an inner sub-Neptune and an outer warm Saturn [Mantovan et al. 2024a](https://ui.adsabs.harvard.edu/abs/2024A%26A...682A.129M/abstract).
- Confronting compositional confusion through the characterization of the sub-Neptune orbiting HD 77946 [Palethorpe et al. 2024](https://ui.adsabs.harvard.edu/abs/2024MNRAS.529.3323P/abstract)
- TASTE. V. A new ground-based investigation of orbital decay in the ultra-hot Jupiter WASP-12b [Leonardi et al. 2024](https://ui.adsabs.harvard.edu/abs/2024A%26A...686A..84L/abstract)
- A new dynamical modeling of the WASP-47 system with CHEOPS observations [Nascimbeni et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...673A..42N/abstract)
- GAPS XXXVII. A precise density measurement of the young ultra-short period planet TOI-1807 b [Nardiello et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract)
- Investigating the architecture and internal structure of the TOI-561 system planets with CHEOPS, HARPS-N, and TESS [Lacedelli et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.511.4551L/abstract)
- Independent validation of the temperate Super-Earth HD79211 b using HARPS-N [DiTommaso et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv221012211D/abstract)
- Multi-mask least-squares deconvolution: extracting RVs using tailored masks [Lienhard et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5328L/abstract)
- Kepler-102: Masses and Compositions for a Super-Earth and Sub-Neptune Orbiting an Active Star [Brinkman et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022arXiv221105196B/abstract)
- TIC 257060897b: An inflated, low-density, hot-Jupiter transiting a rapidly evolving subgiant star [Montalto et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.2908M/abstract)
- A PSF-based Approach to TESS High quality data Of Stellar clusters (PATHOS) - IV. Candidate exoplanets around stars in open clusters: frequency and age-planetary radius distribution [Nardiello et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.3767N/abstract)
- An unusually low density ultra-short period super-Earth and three mini-Neptunes around the old star TOI-561 [Lacedelli et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021MNRAS.501.4148L/abstract)
- K2-111: an old system with two planets in near-resonance [Mortier et al. (2020)](https://ui.adsabs.harvard.edu/abs/2020MNRAS.499.5004M/abstract)
- GAPS XXIX. No detection of reflected light from 51 Peg b using optical high-resolution spectroscopy [Scandariato et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...646A.159S/abstract)
- GAPS XXVIII. A pair of hot-Neptunes orbiting the young star TOI-942 [Carleo et al. (2021)](https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..71C/abstract)

And many more I may have missed!
