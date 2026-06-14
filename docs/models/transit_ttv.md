# Transit models for TTV measurements

The TTV transit models fit independent mid-transit times instead of forcing all
events to follow the same `Tc + N * P` ephemeris.  They keep the usual transit
shape parameters (`R_Rs`, impact parameter or inclination, stellar density or
`a_Rs`, eccentricity parametrization and limb darkening), but replace the single
global conjunction epoch with one fitted `Tc` per observed transit.

These models are the right choice when the measured transit times themselves are
the scientific product.  For ordinary transit fits with one shared
inferior-conjunction epoch, use
[Transit light-curve models](transit_lightcurves.md).

The examples in the external
[PyORBIT_examples TTV directory](https://github.com/LucaMalavolta/PyORBIT_examples/tree/main/ttv_measurement)
show the same four workflows described below: one file per transit, transit
identifiers in a `subset` column, an external transit-time list, and ancillary
columns carrying the transit identifiers.

## Shared planet setup

TTV models force the planet to use inferior-conjunction times.  It is still good
practice to state this explicitly in the planet common object:

```yaml
planet_b:
  common: planets
  orbit: circular
  parametrization: Eastman2013
  use_time_inferior_conjunction: True
```

The period `P` is still part of the planet common object.  It provides the
linear ephemeris used to label or initialize transit events, while the observed
mid-times are fitted as independent `Tc` parameters.

The optional model keywords `Tc_boundaries`, `Tc_bounds` and `Tc_priors` can be
used to override the default bounds or priors assigned to the fitted transit
times.

## Model families

| Workflow | Batman model | PyTransit model | Transit-time parameters |
| --- | --- | --- | --- |
| One dataset per transit | `batman_transit_ttv` | `pytransit_transit_ttv` | one dataset-level `Tc` for each dataset |
| Transit IDs from `subset` | `batman_transit_ttv_subset` | `pytransit_transit_ttv_subset` | `Tc_0`, `Tc_1`, ... from the dataset subsets |
| External transit-time list | `batman_transit_ttv_tclist` | `pytransit_transit_ttv_tclist` | `Tc_N` for events selected from the list |
| Transit IDs from ancillary columns | `batman_transit_ttv_ancillary` | `pytransit_transit_ttv_ancillary` | `Tc_N` from the selected ancillary flag |

The aliases `subset_batman_transit_ttv`, `subset_pytransit_transit_ttv`,
`tclist_batman_transit_ttv`, `tclist_pytransit_transit_ttv`,
`ancillary_batman_transit_ttv` and `ancillary_pytransit_transit_ttv` are also
accepted.

## Sharing transit times across datasets

By default, fitted TTV parameters can be local to each dataset/model instance.
Set `use_shared_ttv: True` when the same transit time must be shared by multiple
light curves, instruments or passbands.

```yaml
lc_model:
  model: batman_transit_ttv_tclist
  planets: [b]
  limb_darkening: ld_quadratic
  use_shared_ttv: True
```

This is especially useful for simultaneous multi-instrument observations and for
multi-planet configurations where each planet has its own transit-time list.

## One dataset per transit

Use `batman_transit_ttv` or `pytransit_transit_ttv` when each input file contains
one transit event.  This is the most explicit layout: every dataset gets its own
dataset-level `Tc`.

```yaml
input:
  LCdata_transit00:
    file: individual_transits/lc_transit00.dat
    kind: Phot
    models:
      - lc_model

  LCdata_transit01:
    file: individual_transits/lc_transit01.dat
    kind: Phot
    models:
      - lc_model

lc_model:
  model: batman_transit_ttv
  planets: [b]
  limb_darkening: ld_quadratic
```

Use this form when the light curves are already split by transit and no extra
transit identifier column is needed.

## TTVs from dataset subsets

Use the `_subset` models when one file contains many transits and the dataset
has a `subset` column identifying the event.  The example data in
`PyORBIT_examples/ttv_measurement/subset_transits` use a header like:

```text
# time flux flux_error jitter offset subset
```

Each non-negative subset value identifies a fitted transit time.  For example,
subset `0` maps to `Tc_0`, subset `1` maps to `Tc_1`, and so on.

```yaml
input:
  LCdata_inst0:
    file: subset_transits/lc_inst0.dat
    kind: Phot
    models:
      - lc_model_inst0

lc_model_inst0:
  model: pytransit_transit_ttv_subset
  planets: [b]
  limb_darkening: ld_quadratic
  use_shared_ttv: True
```

This layout is convenient when different instruments observe the same sequence
of transits and the fitted `Tc_N` values must be common to all of them.

## TTVs from an external transit-time list

Use the `_tclist` models when the expected transit windows are stored in a
separate file.  The examples use a table with a planet label, transit identifier,
expected transit time and transit window:

```text
# planet transit_id transit_time transit_window
b 0 2.170000 0.50000
b 1 6.330000 0.50000
```

Declare the list in the planet common object:

```yaml
planet_b:
  common: planets
  orbit: circular
  parametrization: Eastman2013
  use_time_inferior_conjunction: True
  TTV_Tc_list: tclist_transits/WASP47_simulated_planetb_tclist.dat

lc_model:
  model: batman_transit_ttv_tclist
  planets: [b]
  limb_darkening: ld_quadratic
  minimum_number_of_observations: 20
  use_shared_ttv: False
```

The model searches the dataset around each listed transit window and creates a
fitted `Tc_N` parameter for the events with enough data points.  Increase
`minimum_number_of_observations` to skip poorly covered windows.

For a multi-planet TTV fit, give each planet its own `TTV_Tc_list` and include
all planets in the same model:

```yaml
lc_model:
  model: batman_transit_ttv_tclist
  planets: [b, c]
  limb_darkening: ld_quadratic
  use_shared_ttv: True
```

## TTVs from ancillary transit flags

Use the `_ancillary` models when the input file already includes one column that
labels the transit number for each planet.  The two-planet examples use columns
like:

```text
# time flux flux_error jitter offset subset id_transit_b id_transit_c
```

Connect the planet to the correct ancillary column with `TTV_Tc_flag`:

```yaml
planet_b:
  common: planets
  orbit: circular
  parametrization: Eastman2013
  use_time_inferior_conjunction: True
  TTV_Tc_flag: id_transit_b

lc_model:
  model: pytransit_transit_ttv_ancillary
  planets: [b]
  limb_darkening: ld_quadratic
  use_shared_ttv: False
```

Rows with a negative flag are ignored; non-negative identifiers are mapped to
the corresponding `Tc_N` parameter.

## Choosing the workflow

Use the one-dataset-per-transit models for small, already segmented light-curve
sets.  Use `_subset` when a standard PyORBIT `subset` column already identifies
the transit event.  Use `_tclist` when the event windows are easier to maintain
outside the light-curve files, especially for multi-planet systems.  Use
`_ancillary` when the transit IDs are already part of the data table or when
different planets require different event labels in the same file.
