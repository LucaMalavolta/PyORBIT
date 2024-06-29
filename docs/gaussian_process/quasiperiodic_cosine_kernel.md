(quasiperiodic_cosine_kernel)=

# Quasi-periodic with cosine kernel


This kernel has been introduced by [Perger et al. 2021](https://ui.adsabs.harvard.edu/abs/2021A%26A...645A..58P/abstract). The kernel has been implemented in `PyORBIT` without relying on any other package.

If we define $\tau = t_i-t_j$ :

```{math}
:label: quasiperiodic_cosine

G(\tau ) = \exp{\frac{-2 \tau^2}{P_\mathrm{dec}^2}} *  \left [ H_\mathrm{amp}^2 \exp{-\frac {
  \sin^2{( \pi \tau / P_\mathrm{rot} )}}{2 O_\mathrm{amp} ^2}} + C_\mathrm{amp}^2  \cos \frac{ 4\pi \tau}{P_\mathrm{rot}} \right ]

```

where $P_\mathrm{rot}$  is equivalent to the rotation period of the star, $O_\mathrm{amp}$ is the coherence scale, and $P_\mathrm{dec} $ is usually associated with the decay time scale of the active regions. Within `PyORBIT`, the amplitude of the quasi-periodic part of the kernel $h_1$ and the amplitude of the cosine part $h_2$ have been labeled as `Hamp` and `Camp` respectively.


```{important}
As for the quasi-periodic kernel, mind the possible presence of a factor 2 in the denominator of the aperiodic variation (i.e., $2 \lambda$ rather than $\lambda$)
```


## Model definition and requirements

**model name**: `tinygp_quasiperiodic_cosine`
- required common object: `activity`
- this model relies on `tinygp`

**model name**: `gp_quasiperiodic_cosine`
- required common object: `activity`
- *direct* implementation relying only on `scipy` and `numpy`


## Model parameters

The following parameters will be inherited from the common model (column *Common?: common*) or a different value will be assigned for each dataset (column *Common?: dataset*)

| Name        | Parameter | Common?  | Definition  | Notes |
| :---        | :-------- | :-------------  | :-----  | :---- |
| Prot      | Rotational period of the star $\theta$ | common | ``activity``     | |
| Pdec      | Decay time scale of active regions $\lambda$ | common | ``activity``     | |
| Oamp | Coherence scale $w$ | common | ``activity`` |   |
| Hamp  | Amplitude of the kernel | dataset | ``activity``     | |
| Camp  | Amplitude of the cosine part of the kernel | dataset | ``activity``     | |


## Keywords

Model-wide keywords, with the default value in boldface.

**hyperparameters_condition**
* accepted values: `True` | **`False`**
* activate the conditions $ \lambda ^ 2 > (3/4 \pi) \theta ^2 w ^ 2 $ (adapted from [Rajpaul 2017](https://ui.adsabs.harvard.edu/abs/2017PhDT.......229R/abstract) and [Rajpaul et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1847R/abstract) to take into account the factor 2 in the denominator of the aperiodic variation) to ensure that the QP function has at least one non-trivial turning point.

**rotation_decay_condition**
* accepted values: `True` | **`False`**
* if activated, it ensures that the decay time scale of the activity regions $\lambda$ is at least twice the rotational period of the star $\theta$

**use_stellar_rotation_period**
* accepted values: `True` | **`False`**
* if activated, the parameter `Prot` from the `activity` *common model* will be replaced by the parameter `rotation_period` from the `star_parameters` *common model*. In this way, a unique parameter can be used by different models, e.g., stellar activity and Rossiter-McLaughlin modeling. It can also be useful if you want to use independent GP hyperparameters over several observational seasons while using a single parameter for the rotational period of the star.
