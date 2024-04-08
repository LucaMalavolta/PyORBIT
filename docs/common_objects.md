(common_objects)=

# Common objects

In this section of the configuration file, we include all the physical objects
whose properties are going to be shared among several models, and that are
independent of the datasets under analysis.
For example, the orbital period of planet is a *common parameter* and its value
will be the same regardless of the model selected to analyze the data (e.g. the
planet must have the same period if measured with either transits or radial velocities).
As another example, the rotational period of a star and the decay time scale
of active regions can be considered *common parameters*, while the amplitude
of stellar activity will depend on the dataset under analysis as it changes with wavelength.

```{note}
Not all the physical parameters defined in a *common object* will be employed in the modelling, but only those required by the models listed in the ``model`` section of the configuration file and recalled in the ``model`` subsection of each dataset under the input section
```

```{toctree}
:maxdepth: 1
common_objects/planets
common_objects/star
```
