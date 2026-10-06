(abstract_transit)=

# Abstract Transit Model

Transit parameters and keywords are defined in an abstract class that is then inherited by all the different implementations of the transit model and transit-like modelling (phase curves, Rossiter-McLaughlin).

## Model definition and requirements

The abstract transit class is shared by these data-analysis models:

- ``batman_transit``
- ``pytransit_transit``
- ``batman_transit_secondary_phasecurve``
- ``batman_transit_ttv``
- ``pytransit_transit_ttv``

## Model parameters

The concrete transit model takes orbital and geometric parameters from the `planets` common object and, where applicable, limb-darkening parameters from its common object. The shared integration settings below determine how model fluxes are averaged over finite exposures.

## Keywords

**model**
* accepted values: any of the models listed above
* default: the model with the same name assigned to the section

**planets**:
* list of planets to be modelled, as defined in the ``planets`` common model

**supersample_factor**
Number of points used to sample the average value of the light curve over the entire exposure, in order to take into account the deformation of the light curve due to finite exposure time (see [Kipping 2010](https://ui.adsabs.harvard.edu/abs/2010MNRAS.408.1758K))

**exposure_time**
Exposure time of each data point in the lightcurve, expressed in seconds. The keyword is ignored when ``supersample_factor`` is equal to 1.

**nthreads**
* default is 1, it shouldn't be modified
* number of CPU used by ``batman`` in multiprocessing mode. It seems to cause a clash with ``emcee``, while the better performance of multiprocessing in nested sampling algorithms does not really require the usi of multiple cpu in the computation of the lightcurve model.

## Examples

```yaml
models:
  lc_model:
    kind: batman_transit
    planets: b
    supersample_factor: 1
    exposure_time: 120.00
```
