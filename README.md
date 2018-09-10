# PyORBIT

## `PyORBIT` version 6.0 by Luca Malavolta - September 2018   

**New** Version 6 correctly implements model selection with nested sampling performed by `MultiNest`/`PolyChord`. Because of the way the software has been structured, each model must be contained in a different `yaml` file, i.e. the number of planets is not an hyperparameter but it must be fixed a priori. However, the same `yaml` file can be used to run `emcee`, `MultiNest` or `PolyChord` independently.

**Documentation** Some incomplete documentation is available [here](http://pyorbit.readthedocs.io/). For any doubt, feel free to contact me at luca.malavolta_at_unipd.it, I'll be happy to work out together any problem may arise during installation or usage of this software.

`PyORBIT` handles several kinds of datasets, such as radial velocity (RV), activity indexes, and photometry, to simultaneously characterize the orbital parameters of exoplanets and the noise induced by the activity of the host star. RV computation is performed using either non-interacting Kepler orbits or n-body integration. Stellar activity can be modeled either with sinusoids at the rotational period and its harmonics or gaussian process. Offsets and systematics in measurements from several instruments can be modeled as well. Thanks to the modular approach, new methods for stellar activity modeling or parameter estimation can be easily incorporated into the code.

Most of the information can be found in [Malavolta et al. (2016)](https://ui.adsabs.harvard.edu//#abs/2016A&A...588A.118M/abstract) and [Malavolta et al. (2018)](https://ui.adsabs.harvard.edu//#abs/2018AJ....155..107M/abstract).  
