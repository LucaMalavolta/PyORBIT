(abstract_transit)=

# Abstract Transit Model

Transit parameters and keywords are defined in an abstract class that is then inherited by all the different implementations of the transit model and transit-like modelling (phase curves, Rossiter-McLaughlin)..

Models relying on the Transit abstract model:
- ``batman_transit``
- ``pytransit_batman``
- ``batman_transit_secondary_phasecurve``
- ``batman_transit_with_ttv``
- ``pytransit_transit_with_ttv``


### Keywords

**kind**
* accepted values: any of the models listed above
* default: the model with the same name assigned to the section

**planets**:


**nthreads**
* default is 1, it shouldn't be modified
* number of CPU used by batman in multiprocessing mode. It seems to cause a clash with ``emcee``, while the better performance of multiprocessing in nested sampling algorithms do not really require the usi of multiple cpu in the computation of the lightcurve model.

### Example

```yaml
models:
  lc_model:
    kind: batman_transit
    planets: b
    supersample_factor: 1
    exposure_time: 120.00
```
