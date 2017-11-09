# PyORBIT

**`PyORBIT` v5.0 by Luca Malavolta - 2017**    

*Full documentation coming soon*

PyORBIT handles several kinds of datasets, such as radial velocity (RV), activity indexes, and photometry, to simultaneously characterize the orbital parameters of exoplanets and the noise induced by the activity of the host star. RV computation is performed using either non-interacting Kepler orbits or n-body integration. Stellar activity can be modeled either with sinusoids at the rotational period and its harmonics or gaussian process. Offsets and systematics in measurements from several instruments can be modeled as well. Thanks to the modular approach, new methods for stellar activity modeling or parameter estimation can be easily incorporated into the code.

Most of the information can be found in the paper by [Malavolta et al. (2016)][Malavolta2016]
