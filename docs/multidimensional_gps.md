# Multidimensional GPs

In the original Gaussian process framework ([Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract), [Rajpaul et al. 2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.507.1847R/abstract)) the radial velocity datasets ($\Delta \mathrm{RV}$, after removing the deterministic part) and two activity indicators (in this example, $\mathrm{BIS}$ and  $\log{R^{\prime}_\mathrm{HK}}$) are modeled as a liner combination of an underlying Gaussian process $G(t)$  and its first derivative $G^\prime (t)$.

```{math}
:label: gp_framework_original

\Delta \mathrm{RV} & = V_c G(t) + V_r G^\prime (t) \\
\mathrm{BIS} & = B_c G(t) + B_r G^\prime (t) \\
\log{R^{\prime}_\mathrm{HK}} & = L_c G(t) \\
```

`PyORBIT` implementation produces results that are perfectly consistent with the GP framework by [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract)  and the multidimensional GP by [Barragán et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..866B/abstract), as shown in [Nardiello et al. 2022, appendix D](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract). Note that the signs of the coefficients may vary according to the employed definition of  $\Delta t = (t_i-t_j)$.

```{admonition} Give credits to the authors
If you use the multidimensional GP, remeber to cite [Rajpaul et al. 2015](https://ui.adsabs.harvard.edu/abs/2015MNRAS.452.2269R/abstract)  and [Barragán et al. 2022](https://ui.adsabs.harvard.edu/abs/2022MNRAS.509..866B/abstract).
For more details on the `PyORBIT` implementation of multidimensional GP, please refer to [Nardiello et al. 2022, appendix D](https://ui.adsabs.harvard.edu/abs/2022A%26A...664A.163N/abstract)

```

```{toctree}
:maxdepth: 1
multidimensional_gps/md_quasiperiodic.md
multidimensional_gps/md_spleaf.md

```
