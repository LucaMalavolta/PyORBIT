(gp_rv_lc_fit)=

# Simultaneous Gaussian Process regression on light curves and radial velocities

Imagine we have both photometric and spectroscopic time series of a rather active star. Photometry will provide essential infomration on the rotatinoal period of the star and the decay time scale of activity regions, while spectroscopic activity indexes are essential to disentagle planetary and stellar signal in radial velocities. 

In this example, we will analyze the planetary system around K2-141 [Malavolta et al. 2018](https://ui.adsabs.harvard.edu/abs/2018AJ....155..107M/abstract). The spectroscopic dataset is taken from [Bonomo et al. 2023](https://ui.adsabs.harvard.edu/abs/2023A%26A...677A..33B/abstract) and includes Radial velocity data, Bisector Inverse Span, and S-index, all gathered with HARPS-N. 

Photometric data comes from different sources:
- K2 campaign 12, corrected fro systematics as described in [Vanderburg and Johnson 2014](https://ui.adsabs.harvard.edu/abs/2014PASP..126..948V/abstract) and available at the [MAST K2SFF website](https://archive.stsci.edu/hlsp/k2sff)
- TESS sector 42, [info available here](https://tess.mit.edu/observations/sector-42/)
- TESS sector 70, [info available here](https://tess.mit.edu/observations/sector-70/)

```{figure} plots/K2-141_lighcurves.png
:alt: lightcurve of K2-141 
:width: 80 %
Photometric data for K2-141. The stellar activity modulation is clearly identifiable.
```

We want to model stellar activity in photometry simultaneously with the spectrosocpic time series. 

- For photometric data, we use the [exponential-sine periodic kernel] (ESP) implemented in 'S+LEAF'. Each dataset will have its independent covariance matrix.
- for spectroscopic data, we use the [multidimensional quasi-periodic kernel](../multidimensional_gps/md_quasiperiodic.md)
 through `tinygp`.

For the photometric data, we also include:
- a normalization factor, independent for each dataset
- the transit model for the two planets. 

For K2 and TESS we use two different transit models, as they have different limb darkening coefficients. The `lc_model_k2` and `lc_model_tess` models differs in the limb darkening coefficients, `limb_darkening_k2` and `limb_darkening_tess`, specififed in the `common:star` section and recalled in `models:lc_model_k2:limb_darkening` and models:lc_model_tess:limb_darkening`, respectively. 

For the spectroscopic data, we only have to include the `radial_velocities` model due to the planets in the radial velocity dataset. 

The input section takes this form:


```{code-block} yaml
:lineno-start: 1
inputs:
  LC_K2_c12:
    file: datasets/K2-141_K2_KSFF_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_k2
      - spleaf_esp
      - normalization_factor
  LC_TESS_s42:
    file:  datasets/K2-141_TESS_s42_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_tess
      - spleaf_esp
      - normalization_factor
  LC_TESS_s70:
    file:  datasets/K2-141_TESS_s70_PyORBIT.dat
    kind: Phot
    models:
      - lc_model_tess
      - spleaf_esp
      - normalization_factor
  RVdata:
    file: datasets/K2-141_RV_PyORBIT.dat
    kind: RV
    models:
      - radial_velocities
      - gp_multidimensional
  BISdata:
    file: datasets/K2-141_BIS_PyORBIT.dat
    kind: BIS
    models:
      - gp_multidimensional
  Sdata:
    file: datasets/K2-141_Sindex_PyORBIT.dat
    kind: S_index
    models:
      - gp_multidimensional
```


