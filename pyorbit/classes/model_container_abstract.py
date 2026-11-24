from pyorbit.subroutines.common import *
from packaging.version import parse as parse_version
import pyorbit
import time

__all__ = ["ModelContainer"]

# F-batch (Obs d): the ``lnprob_batch`` per-theta fallback can be hit on
# EVERY batched call for a mis-gated container; log the first occurrence and
# every 100th, stay quiet otherwise (keeps the failure path diagnosable
# without flooding stdout).
_LNPROB_BATCH_FALLBACK_COUNT = {"n": 0}
_LNPROB_BATCH_FALLBACK_PERIOD = 100


class ModelContainer(object):
    def __init__(self):

        self.pyorbit_version = pyorbit.__version__

        # FIXME: change the renaming and link directly the variables in constants module

        # These dictionaries are dedicated to RVs and T0s computed via dynamical integration
        # Since the individual contribution of each planet cannot be disentagled, the output of
        # the dynamical integration cannot be associated to a model.
        self.dynamical_dict = {}
        self.dynamical_model = None

        self.dataset_dict = {}

        # Dictionaries including the common models and data-specific models
        self.common_models = {}
        self.models = {}
        # TODO: Added in PyORBIT version 12 beta
        self.parent_models = {}

        self.include_priors = True

        self.ndata = None
        self.ndof = None

        self.starting_point = None
        self.starting_point_flag = False
        self.recenter_bounds_flag = True

        self.bounds = None
        self.spaces = None
        self.priors = None
        self.range = None
        self.ndim = 0
        self.pam_names = ""

        self.Tref = 2460000.000

        self.ordered_planets = {}

    def model_setup(self):

        # GPU gate (Round 1): fail fast, before any model is initialized, if the
        # fit contains a model that is not JAX-native. This is a strict no-op
        # while the GPU is disabled (the default), so the CPU path is untouched.
        # The check is by *class* (see pyorbit.common.gpu_backend.gate), which is
        # robust to PyORBIT's internal model naming; under `gpu: true` a rejected
        # model raises GPUSupportedModelError naming the JAX equivalent to use.
        pyorbit.common.gpu_backend.gate(self.models, self.common_models)

        # Fixing back-compatibility issues within  v8 sub-versions
        if not hasattr(self, "ordered_planets"):
            self.ordered_planets = {}

        # First step: setting up the correct associations between models and dataset

        for model_name, model in self.common_models.items():
            try:
                model_conf = model.model_conf
                if model_conf is None:
                    model_conf = {}
            except:
                model_conf = {}

            model.initialize_model(self, **model_conf)

        for model_name, model in self.models.items():
            try:
                model_conf = model.model_conf
                if model_conf is None:
                    model_conf = {}
            except:
                model_conf = {}

            model.initialize_model(self, **model_conf)

            # TODO: added in PyORBIT version 12 beta
            model.initialize_model_parameters(self, **model_conf)

            model.change_parameter_status(self, **model_conf)

            for dataset_name in list(
                OrderedSet(model_conf) & OrderedSet(self.dataset_dict)
            ):
                model.initialize_model_dataset(
                    self, self.dataset_dict[dataset_name], **model_conf
                )

        # TODO: added in PyORBIT version 12 beta
        for model_name, model in self.parent_models.items():
            try:
                model_conf = model.model_conf
                if model_conf is None:
                    model_conf = {}
            except:
                model_conf = {}
            model.initialize_model(self, **model_conf)

            for dataset_name in list(
                OrderedSet(model_conf) & OrderedSet(self.dataset_dict)
            ):
                model.initialize_model_dataset(
                    self, self.dataset_dict[dataset_name], **model_conf
                )

        if self.dynamical_model:
            self.dynamical_model.to_be_initialized = True

    def print_model_info(self):
        for model_name, model in self.common_models.items():
            if model.skip_print_info:
                continue
            model.print_info()

        for model_name, model in self.models.items():
            if model.skip_print_info:
                continue
            model.print_info()

        for model_name, model in self.parent_models.items():
            if model.skip_print_info:
                continue
            model.print_info()

    def boundaries_setup(self):
        # This routine setup the boundary array and at the same time
        # define a dictionary with the name of the arrays and their
        # positions in bounds/theta array so that they can be accessed
        # without using nested counters

        self.ndim = 0
        output_lists = {
            "bounds": [],
            "spaces": [],
            "priors": [],
        }

        for model_name, model in self.models.items():
            if len(model.common_ref) > 0:
                for common_ref in model.common_ref:
                    model.default_bounds.update(
                        self.common_models[common_ref].default_bounds
                    )
                    model.default_spaces.update(
                        self.common_models[common_ref].default_spaces
                    )
                    model.default_priors.update(
                        self.common_models[common_ref].default_priors
                    )
                    self.ndim, output_lists = self.common_models[
                        common_ref
                    ].define_parameter_properties(
                        self.ndim, output_lists, model.list_pams_common
                    )
            else:
                pass

        for dataset_name, dataset in self.dataset_dict.items():
            self.ndim, output_lists = dataset.define_parameter_properties(
                self.ndim, output_lists, dataset.list_pams
            )

            for model_name in dataset.models:
                self.ndim, output_lists = self.models[
                    model_name
                ].define_parameter_properties(self.ndim, output_lists, dataset.name_ref)

        self.bounds = np.asarray(output_lists["bounds"])
        self.spaces = output_lists["spaces"]
        self.priors = output_lists["priors"]
        self.range = self.bounds[:, 1] - self.bounds[:, 0]

        if any(self.range < 0.0):
            print("UNRECOVERABLE ERROR from configuration file")
            print("    Some parameter has lower bound greater than upper bound")
            print()
            quit()

    def initialize_logchi2(self):

        # Second step: define the number of variables and setting up the boundaries
        # To be done only the first time ever
        self.ndata = 0
        for dataset_name, dataset in self.dataset_dict.items():
            if not dataset.models:
                continue
            self.ndata += dataset.n
        self.ndof = self.ndata - self.ndim

    def starting_points_setup(self):

        self.starting_point = np.average(self.bounds, axis=1)
        self.starting_point[:] = None

        for model_name, model in self.common_models.items():
            model.define_starting_point(self.starting_point)

        for dataset_name, dataset in self.dataset_dict.items():
            dataset.define_starting_point(self.starting_point)

            for model in dataset.models:
                self.models[model].define_starting_point(
                    self.starting_point, dataset_name
                )

    def check_bounds(self, theta):

        for ii in range(0, self.ndim):
            if not (self.bounds[ii, 0] < theta[ii] < self.bounds[ii, 1]):
                return False

        period_storage = []
        period_storage_ordered = np.zeros(len(self.ordered_planets))

        for model_name, model in self.common_models.items():
            if model.model_class == "planet":
                """ Step 1: retrieve the planet period"""
                period = model.transformation["P"](
                    theta, model.fixed, model.parameter_index["P"]
                )

                """ Step 2: save the all planet periods into a list"""
                period_storage.extend([period])

                """ Step 3: save the period of the planet in the ordered list"""
                if model_name in self.ordered_planets:
                    period_storage_ordered[self.ordered_planets[model_name]] = period
                    # print('    ', model_name, self.ordered_planets[model_name], period_storage_ordered)

                """ Step 4: check if the eccentricity is within the given range"""
                if "e" in model.parameter_index:
                    e = model.transformation["e"](
                        theta, model.fixed, model.parameter_index["e"]
                    )

                    if not model.bounds["e"][0] <= e < model.bounds["e"][1]:
                        # print('eccentricity>1')
                        # print()
                        return False

                """ Step 5: check if the impact parameter is below 1 + Rp/Rs """
                if "b" in model.parameter_index and "R_Rs" in model.parameter_index:
                    b = model.transformation["b"](
                        theta, model.fixed, model.parameter_index["b"]
                    )
                    R = model.transformation["R_Rs"](
                        theta, model.fixed, model.parameter_index["R_Rs"]
                    )
                    if not np.abs(b) <= 1 + R:
                        return False

                """ Step 6 eclipse depth must be greater than the amplitude of
                phase variation (during the eclipse we only have the stellar line) """
                if (
                    "phase_amp" in model.parameter_index
                    and "delta_occ" in model.parameter_index
                ):
                    phase_amp = model.transformation["delta_occ"](
                        theta, model.fixed, model.parameter_index["phase_amp"]
                    )
                    delta_occ = model.transformation["delta_occ"](
                        theta, model.fixed, model.parameter_index["delta_occ"]
                    )
                    if phase_amp > delta_occ:
                        return False

        # """ Step 6 eclipse depth must be greater than the amplitude of
        # phase variation (during the eclipse we only have the stellar line) """
        # for dataset_name, dataset in self.dataset_dict.items():
        #    if 'phase_amp' in dataset.parameter_index and 'delta_occ' in dataset.parameter_index:
        #        phase_amp = dataset.transformation['delta_occ'](theta,
        #                                        dataset.fixed,
        #                                        dataset.parameter_index['phase_amp'])
        #        delta_occ = dataset.transformation['delta_occ'](theta,
        #                                        dataset.fixed,
        #                                        dataset.parameter_index['delta_occ'])
        #        if phase_amp > delta_occ: return False

        """ Step 7 check for overlapping periods (within 2.5% arbitrarily chosen)"""
        for i_n, i_v in enumerate(period_storage):
            if i_n == len(period_storage) - 1:
                break
            if np.amin(np.abs(period_storage[i_n + 1 :] - i_v)) / i_v < 0.025:
                # print('overlapping periods  detected')
                # print()
                return False

        """ Step 8 check if the planet are ordered"""
        for i_n, i_v in enumerate(period_storage_ordered):
            if i_n == len(period_storage_ordered) - 1:
                break
            if np.amin(period_storage_ordered[i_n + 1 :] - i_v) < 0.0:
                # print('inverted order detected')
                # print()
                return False

        return True

    def __call__(self, theta, include_priors=True):

        log_priors, log_likelihood = self.log_priors_likelihood(theta)

        if self.include_priors and include_priors:
            return log_priors + log_likelihood
        else:
            return log_likelihood

    def log_priors_likelihood(self, theta, return_priors=True):

        log_priors = 0.00
        log_likelihood = 0.00

        """
        Constant term added either by dataset.model_logchi2() or gp.log_likelihood()
        """

        if not self.check_bounds(theta):
            if return_priors is False:
                return -np.inf
            else:
                return -np.inf, -np.inf

        if self.dynamical_model is not None:
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dynamical_output = self.dynamical_model.compute(self, theta)

            # print(dynamical_output['stable'], dynamical_output['pass'], dynamical_output['tc_check'])

            if not (dynamical_output["stable"] and dynamical_output["pass"]):
                if return_priors is False:
                    return -np.inf
                else:
                    return -np.inf, -np.inf

            if not dynamical_output["tc_check"]:
                if return_priors is False:
                    return -np.inf
                else:
                    return -np.inf, -np.inf

        for model_name, model in self.common_models.items():
            log_priors += model.return_priors(theta)

        delayed_lnlk_computation = []
        residuals_analysis = {}

        for dataset_name, dataset in self.dataset_dict.items():
            logchi2_gp_model = None
            compute_gp_residuals = False

            dataset.model_reset()
            parameter_values = dataset.convert(theta)
            dataset.compute(parameter_values)

            log_priors += dataset.return_priors(theta)

            if "none" in dataset.models or "None" in dataset.models:
                continue
            if not dataset.models:
                continue

            # TODO: Added in PyORBIT version 12 beta
            multiple_planets_models = {}
            skip_loglikelihood = False

            for model_name in dataset.models:
                log_priors += self.models[model_name].return_priors(theta, dataset_name)

                parameter_values = {}

                for common_ref in self.models[model_name].common_ref:
                    parameter_values.update(
                        self.common_models[common_ref].convert(theta)
                    )

                # TODO: here, possible check for models requiring a secondary star
                if getattr(self.models[model_name], "requires_secondary_star", False):
                    binary_ref = self.models[model_name].binary_ref
                    parameter_values.update(
                        self.common_models[binary_ref].convert_with_name(theta, "bin")
                    )

                parameter_values.update(
                    self.models[model_name].convert(theta, dataset_name)
                )

                # ** NEW in PyORBIT version 12
                # dataset independent calculation for the model
                # The dataset is passed only for keywords such as reference time and so

                # TODO: Added in PyORBIT version 12 beta
                if getattr(self.models[model_name], "accept_multiple_planets", False):
                    parent_model = self.models[model_name].parent_model
                    planet_ref = self.models[model_name].planet_ref
                    parameter_values.update(
                        self.common_models[planet_ref].convert_with_name(
                            theta, planet_ref
                        )
                    )
                    self.models[model_name].update_parameter_values(
                        parameter_values, planet_ref
                    )

                    self.parent_models[parent_model].parameter_values.update(
                        parameter_values
                    )

                    try:
                        multiple_planets_models[parent_model].append(planet_ref)
                    except KeyError:
                        multiple_planets_models[parent_model] = [planet_ref]
                    continue
                else:
                    self.models[model_name].precompute(parameter_values, dataset)

                """ compute the log-likelihood for those models using external datasets,
                    but only when multiple planets are not involved, otherwise the log-likelihood is computed at the parent model level"""
                if getattr(self.models[model_name], "external_dataset", False):
                    skip_loglikelihood = True
                    log_likelihood += self.models[model_name].compute_loglikelihood(
                        parameter_values, dataset
                    )
                    continue

                """ residuals will be computed following the definition in Dataset class
                    This section has never been tested
                """
                if getattr(self.models[model_name], "residuals_analysis", False):
                    skip_loglikelihood = True
                    if self.models[model_name].gp_before_correlation:
                        compute_gp_residuals = True

                    try:
                        residuals_analysis[model_name]["parameter_values"] = (
                            parameter_values
                        )
                    except:
                        residuals_analysis[model_name] = {
                            "parameter_values": parameter_values
                        }

                    if dataset_name == self.models[model_name].x_dataset:
                        residuals_dataset_label = "x"
                        residuals_analysis[model_name]["x_gp_model"] = None
                        residuals_analysis[model_name]["x_gp_parameters"] = None
                    else:
                        # Subtract the priors previously added, to avoid including the priors twice
                        log_priors -= self.models[model_name].return_priors(
                            theta, dataset_name
                        )

                        residuals_dataset_label = "y"
                        residuals_analysis[model_name]["y_gp_model"] = None
                        residuals_analysis[model_name]["y_gp_parameters"] = None
                    continue

                if getattr(self.models[model_name], "internal_likelihood", False):
                    logchi2_gp_model = model_name
                    continue

                # if getattr(self.models[model_name], 'model_class', None) is 'common_jitter':
                if getattr(self.models[model_name], "jitter_model", False):
                    dataset.jitter += self.models[model_name].compute(
                        parameter_values, dataset
                    )
                    continue

                # print('    Computing model: ', model_name, ' for dataset: ', dataset_name)
                # print(dataset_name, dataset.kind, dataset.dynamical, self.models[model_name].dynamical_model)

                if getattr(dataset, "dynamical", False) and getattr(
                    self.models[model_name], "dynamical_model", False
                ):
                    if dataset.kind == "photometry":
                        dataset.buffer_model = self.models[
                            model_name
                        ].compute_dynamical(
                            parameter_values, dataset, dynamical_output[dataset_name]
                        )
                    else:
                        dataset.external_model = dynamical_output[dataset_name]

                if dataset.normalization_model is None and (
                    self.models[model_name].unitary_model
                    or self.models[model_name].normalization_model
                ):
                    dataset.normalization_model = np.ones(dataset.n, dtype=np.double)

                if self.models[model_name].unitary_model:
                    # print('CHOICE D ', model_name, self.models[model_name].random_number, np.shape(dataset.buffer_model), np.shape(dataset.unitary_model))
                    dataset.unitary_model += self.models[model_name].compute(
                        parameter_values, dataset
                    )
                elif self.models[model_name].normalization_model:
                    dataset.normalization_model *= self.models[model_name].compute(
                        parameter_values, dataset
                    )
                else:
                    dataset.additive_model += self.models[model_name].compute(
                        parameter_values, dataset
                    )

            # TODO: Added in PyORBIT version 12 beta
            for parent_model in multiple_planets_models:
                planet_list = multiple_planets_models[parent_model]

                self.parent_models[parent_model].precompute(planet_list, dataset)

                if getattr(self.parent_models[parent_model], "external_dataset", False):
                    skip_loglikelihood = True
                    log_likelihood += self.parent_models[
                        parent_model
                    ].compute_loglikelihood(planet_list, dataset)
                    continue

                model_out = self.parent_models[parent_model].compute(
                    planet_list, dataset
                )

                if dataset.normalization_model is None and (
                    self.parent_models[parent_model].unitary_model
                    or self.parent_models[parent_model].normalization_model
                ):
                    dataset.normalization_model = np.ones(dataset.n, dtype=np.double)

                if self.parent_models[parent_model].unitary_model:
                    # print('CHOICE D ', model_name, self.models[model_name].random_number, np.shape(dataset.buffer_model), np.shape(dataset.unitary_model))
                    dataset.unitary_model += model_out
                elif self.parent_models[parent_model].normalization_model:
                    dataset.normalization_model *= model_out
                else:
                    dataset.additive_model += model_out

            dataset.compute_model()
            dataset.compute_residuals()

            """ Gaussian Process check MUST be the last one or the program will fail
            that's because for the GP to work we need to know the _deterministic_ part of the model
            (i.e. the theoretical values you get when you feed your model with the parameter values) """

            if logchi2_gp_model:
                parameter_values = {}
                for common_ref in self.models[logchi2_gp_model].common_ref:
                    parameter_values.update(
                        self.common_models[common_ref].convert(theta)
                    )

                parameter_values.update(
                    self.models[logchi2_gp_model].convert(theta, dataset_name)
                )

                """ GP Log-likelihood is not computed now because a single matrix must be
                    computed with the joint dataset"""
                if hasattr(self.models[logchi2_gp_model], "delayed_lnlk_computation"):
                    self.models[logchi2_gp_model].add_internal_dataset(
                        parameter_values, dataset
                    )
                    if logchi2_gp_model not in delayed_lnlk_computation:
                        delayed_lnlk_computation.append(logchi2_gp_model)
                elif skip_loglikelihood:
                    residuals_analysis[model_name][
                        residuals_dataset_label + "_gp_model"
                    ] = logchi2_gp_model
                    residuals_analysis[model_name][
                        residuals_dataset_label + "_gp_parameters"
                    ] = parameter_values
                else:
                    log_likelihood += self.models[logchi2_gp_model].lnlk_compute(
                        parameter_values, dataset
                    )

                if compute_gp_residuals:
                    dataset.residuals_for_regression -= self.models[
                        logchi2_gp_model
                    ].sample_predict(parameter_values, dataset)
                    residuals_analysis[model_name][
                        residuals_dataset_label + "_gp_model"
                    ] = None

            elif not skip_loglikelihood:
                log_likelihood += dataset.model_logchi2()

        """ Correlation between residuals of two datasets through orthogonal distance regression
            it must be coomputed after any other model has been removed from the
            independent variable, if not provided as ancillary dataset
        """
        for model_name in residuals_analysis:
            parameter_values = residuals_analysis[model_name]["parameter_values"]
            x_dataset = self.dataset_dict[self.models[model_name].x_dataset]
            y_dataset = self.dataset_dict[self.models[model_name].y_dataset]
            modelout_xx, modelout_yy = self.models[model_name].compute(
                parameter_values, x_dataset, y_dataset
            )

            x_dataset.residuals -= modelout_xx
            if residuals_analysis[model_name]["x_gp_model"]:
                logchi2_gp_model = residuals_analysis[model_name]["x_gp_model"]
                parameter_values = residuals_analysis[model_name]["x_gp_parameters"]
                log_likelihood += self.models[logchi2_gp_model].lnlk_compute(
                    parameter_values, x_dataset
                )
            else:
                log_likelihood += x_dataset.model_logchi2()

            y_dataset.residuals -= modelout_yy
            if residuals_analysis[model_name]["y_gp_model"]:
                logchi2_gp_model = residuals_analysis[model_name]["y_gp_model"]
                parameter_values = residuals_analysis[model_name]["y_gp_parameters"]
                log_likelihood += self.models[logchi2_gp_model].lnlk_compute(
                    parameter_values, y_dataset
                )
            else:
                log_likelihood += y_dataset.model_logchi2()

        """ In case there is more than one GP model """
        for logchi2_gp_model in delayed_lnlk_computation:
            log_likelihood += self.models[logchi2_gp_model].lnlk_compute()

        """ check for finite log_priors and log_likelihood"""
        if np.isnan(log_priors) or np.isnan(log_likelihood):
            log_likelihood = -np.inf
            log_priors = -np.inf

        for model_name in self.models:
            self.models[model_name].reset_precompute()

        for parent_model in multiple_planets_models:
            self.parent_models[parent_model].reset_precompute()

        # if self.dynamical_model is not None:
        #    if not (dynamical_output['stable'] and dynamical_output['pass']):
        #        log_likelihood = -np.inf
        #        log_priors = -np.inf

        if return_priors is False:
            return log_likelihood
        else:
            return log_priors, log_likelihood

    # ------------------------------------------------------------------
    # F-batch: container-level M-theta numeric lnprob (GPU only)
    # ------------------------------------------------------------------

    def _lnprob_batch_supported(self):
        """F-batch: static structural eligibility for the batched numeric
        lnprob.

        Returns (True, None) only when the container's *numeric* lnprob is
        built exclusively from elements with a proven batch counterpart
        (models exposing the proven ``compute_batch``, dataset-level
        jitter/offset variables). Anything else -- N-body dynamics, GP /
        internal likelihood, external dataset, residuals correlation,
        multi-planet parent models, jitter models, or a model without
        ``compute_batch`` -- is reported by name so the caller can log it
        while falling back to the exact per-theta path. Purely static (no
        theta involved).
        """
        if self.dynamical_model is not None:
            return False, "dynamical_model (N-body)"
        for dataset_name, dataset in self.dataset_dict.items():
            if not dataset.models:
                continue
            if "none" in dataset.models or "None" in dataset.models:
                return False, 'dataset {0}: "none" model entry'.format(dataset_name)
            for model_name in dataset.models:
                model = self.models[model_name]
                if getattr(model, "accept_multiple_planets", False):
                    return False, "{0}: multi-planet/parent model".format(model_name)
                if getattr(model, "external_dataset", False):
                    return False, "{0}: external_dataset".format(model_name)
                if getattr(model, "internal_likelihood", False):
                    return False, "{0}: GP/internal likelihood".format(model_name)
                if getattr(model, "residuals_analysis", False):
                    return False, "{0}: residuals correlation".format(model_name)
                if getattr(model, "jitter_model", False):
                    return False, "{0}: jitter model".format(model_name)
                if not callable(getattr(model, "compute_batch", None)):
                    return False, "{0}: no callable compute_batch".format(model_name)
        return True, None

    def lnprob_batch(self, theta_list):
        """F-batch: M thetas -> (M,) lnprob in one batched GPU evaluation.

        Container-level extension of the proven ``compute_batch`` pattern
        (``pyorbit.models.keplerian_modelling.RVkeplerian.compute_batch``):
        the full *numeric* lnprob (priors + likelihood) of the M thetas is
        evaluated as one batch -- the model curves come from the proven
        single-dispatch (M, N) evaluation and the residuals / chi2
        arithmetic runs on the same (M, N) block on-device -- instead of M
        separate single-theta ``log_priors_likelihood`` calls, each paying
        its own host<->device dispatch round trips.

        Cardinal rule (the batch is a *pure optimization*; it must never
        break a fit):
          * ``gpu: false`` or M <= 1 -> the exact, unchanged per-theta path.
          * out-of-bounds thetas -> -inf, exactly the same thetas as the
            single-theta path (``check_bounds`` is evaluated per theta in a
            plain Python loop; data-dependent guards stay out of the GPU
            body by construction).
          * a container structure without a proven batch counterpart, or
            ANY failure of the batched evaluation (un-vmap-able term, host
            sync inside the body, ...), -> the exact per-theta fallback,
            with a one-line log naming the triggering term. This method
            never raises to its caller because of the batch.
        """
        import pyorbit.common.gpu_backend as gpu_backend

        thetas = [np.asarray(theta, dtype=np.double) for theta in list(theta_list)]
        M = len(thetas)
        if M == 0:
            return np.zeros(0, dtype=np.double)

        def _per_theta(reason):
            """Exact production path (one theta at a time) + fallback log."""
            out = np.zeros(M, dtype=np.double)
            for i, theta in enumerate(thetas):
                if self.include_priors:
                    out[i] = self(theta)
                else:
                    out[i] = self.log_priors_likelihood(theta, return_priors=False)
            _LNPROB_BATCH_FALLBACK_COUNT["n"] += 1
            _n = _LNPROB_BATCH_FALLBACK_COUNT["n"]
            if _n == 1 or _n % _LNPROB_BATCH_FALLBACK_PERIOD == 0:
                _tail = (
                    ""
                    if _n == 1
                    else (
                        "  [repeats throttled: occurrence {0}, next reminder at {1}]".format(
                            _n, _n + _LNPROB_BATCH_FALLBACK_PERIOD
                        )
                    )
                )
                print("lnprob_batch: per-theta fallback ({0}){1}".format(reason, _tail))
            return out

        if not gpu_backend.gpu_enabled():
            return _per_theta("gpu disabled -- per-theta path (gpu: false unchanged)")
        if M == 1:
            return _per_theta("M == 1 -- nothing to batch")

        ok, why = self._lnprob_batch_supported()
        if not ok:
            return _per_theta("not batchable: {0}".format(why))

        """ Step 1 -- domain, per theta, in a plain Python loop (cheap, no
        GPU): the invalid thetas become -inf exactly as in the single-theta
        path. """
        valid = np.asarray([self.check_bounds(theta) for theta in thetas])
        if not valid.any():
            return np.full(M, -np.inf, dtype=np.double)

        """ Step 2 -- numeric lnprob (priors + likelihood) for the valid
        thetas in one batched GPU evaluation; the single host
        materialisation happens when the (Mv,) result is np.asarray-ed
        below (the proven compute_batch pattern). """
        valid_idx = np.nonzero(valid)[0]
        thetas_valid = [thetas[i] for i in valid_idx]

        try:
            priors_v = self._lnprob_batch_priors(thetas_valid)
            likelihood_v = self._lnprob_batch_likelihood(thetas_valid)
        except Exception as exc:
            msg = str(exc)
            return _per_theta(
                "batched evaluation failed: {0}: {1}".format(
                    type(exc).__name__, msg.splitlines()[0] if msg else ""
                )
            )

        """ Same final NaN guard as log_priors_likelihood: NaN -> -inf. """
        total = (
            priors_v if self.include_priors else np.zeros_like(priors_v)
        ) + likelihood_v
        nan = np.isnan(total)
        if nan.any():
            total[nan] = -np.inf

        out = np.full(M, -np.inf, dtype=np.double)
        out[valid_idx] = total
        return out

    def _lnprob_batch_priors(self, thetas):
        """F-batch: log-priors of a list of thetas.

        The *same* ``return_priors`` calls the single-theta path uses
        (common models, dataset, per-dataset models), in the same order and
        applied to each theta in a plain Python loop: closed-form scalar
        arithmetic, no GPU work, so parity holds by construction and there
        is no vmap risk here.
        """
        M = len(thetas)
        prior = np.zeros(M, dtype=np.double)

        for model_name, common_model in self.common_models.items():
            for i, theta in enumerate(thetas):
                prior[i] += common_model.return_priors(theta)

        for dataset_name, dataset in self.dataset_dict.items():
            for i, theta in enumerate(thetas):
                prior[i] += dataset.return_priors(theta)

            if not dataset.models:
                continue
            if "none" in dataset.models or "None" in dataset.models:
                continue

            for i, theta in enumerate(thetas):
                for model_name in dataset.models:
                    prior[i] += self.models[model_name].return_priors(
                        theta, dataset_name
                    )

        return prior

    def _lnprob_batch_likelihood(self, thetas):
        """F-batch: log-likelihood of a list of (in-bounds) thetas.

        Mirrors the single-theta ``log_priors_likelihood`` exactly, lifted
        onto the batch axis:
          * each per-model curve is the proven ``compute_batch`` -- ONE
            jitted (M, N) dispatch, with the per-theta values entering as
            (M,) argument columns, so the cached XLA program is reused for
            any new theta;
          * the per-dataset accumulation (additive / unitary /
            normalization / jitter / offset), residuals and chi2 are the
            same ``dataset.compute`` / ``compute_model`` /
            ``compute_residuals`` / ``model_logchi2`` (GPU branches)
            expressed on the (M, N) block with JAX, in the same order.
        No host sync happens here; the caller materialises once.
        Only called through :meth:`lnprob_batch`, i.e. after the
        structural gate of :meth:`_lnprob_batch_supported`.
        """
        import jax.numpy as jnp

        M = len(thetas)
        loglik = np.zeros(M, dtype=np.double)

        for dataset_name, dataset in self.dataset_dict.items():
            if not dataset.models:
                continue
            if "none" in dataset.models or "None" in dataset.models:
                continue

            n = dataset.n

            """ dataset.convert + dataset.compute (GPU branch): the
            per-theta jitter / offset variable values, applied on the
            (M, N) block in the same order as the single-theta path. """
            pv_list = [dataset.convert(theta) for theta in thetas]

            jitter_block = None
            offset_block = None
            for var in dataset.list_pams:
                generic = dataset.variable_expanded[var]
                if generic == "jitter" and var in dataset.jitter_mask_dev:
                    vals = jnp.asarray([pv[var] for pv in pv_list], dtype=jnp.float64)
                    term = dataset.jitter_mask_dev[var] * vals[:, None]
                    jitter_block = term if jitter_block is None else jitter_block + term
                elif generic == "offset" and var in dataset.offset_mask_dev:
                    vals = jnp.asarray([pv[var] for pv in pv_list], dtype=jnp.float64)
                    term = dataset.offset_mask_dev[var] * vals[:, None]
                    offset_block = term if offset_block is None else offset_block + term

            """ One curve block per model, in dataset order; each is a
            single jitted (M, N) dispatch (the proven pattern). """
            additive_blocks = []
            unitary_blocks = []
            norm_blocks = []
            for model_name in dataset.models:
                model = self.models[model_name]

                parameter_values_list = []
                for theta in thetas:
                    parameter_values = {}
                    for common_ref in model.common_ref:
                        parameter_values.update(
                            self.common_models[common_ref].convert(theta)
                        )
                    if getattr(model, "requires_secondary_star", False):
                        parameter_values.update(
                            self.common_models[model.binary_ref].convert_with_name(
                                theta, "bin"
                            )
                        )
                    parameter_values.update(model.convert(theta, dataset_name))
                    parameter_values_list.append(parameter_values)

                curves = jnp.asarray(
                    model.compute_batch(parameter_values_list, dataset)
                )

                if model.unitary_model:
                    unitary_blocks.append(curves)
                elif model.normalization_model:
                    norm_blocks.append(curves)
                else:
                    additive_blocks.append(curves)

            """ dataset.compute_model (GPU branch) on the (M, N) block, in
            the exact single-theta semantics: a normalization is initialized
            to ones when ANY unitary/normalization model is present, so a
            unitary-only contribution is (1 + sum(u)). """
            additive = None
            for block in additive_blocks:
                additive = block if additive is None else additive + block
            if additive is None:
                additive = jnp.zeros((M, n), dtype=jnp.float64)
            if offset_block is not None:
                additive = additive + offset_block

            unitary = None
            for block in unitary_blocks:
                unitary = block if unitary is None else unitary + block

            norm = None
            for block in norm_blocks:
                norm = block if norm is None else norm * block

            if norm is not None:
                model_m = (
                    additive + (1.0 + (unitary if unitary is not None else 0.0)) * norm
                )
            elif unitary is not None:
                model_m = additive + (1.0 + unitary)
            else:
                model_m = additive

            """ dataset.compute_residuals + dataset.model_logchi2 (GPU
            branch), on the (M, N) block. """
            residuals_m = dataset.y_dev - model_m
            if jitter_block is None:
                env_m = 1.0 / (dataset.e_dev**2.0)
            else:
                env_m = 1.0 / (dataset.e_dev**2.0 + jitter_block**2.0)

            loglik_i = -0.5 * (
                dataset.n * np.log(2.0 * np.pi)
                + jnp.sum(residuals_m**2 * env_m - jnp.log(env_m), axis=1)
            )
            loglik += np.asarray(loglik_i, dtype=np.double)

        return loglik

    def rv_count_ndata_ndim(self):

        ndata = 0
        ndim = 0

        theta_flag = np.zeros(self.ndim, dtype=bool)

        for dataset_name, dataset in self.dataset_dict.items():
            logchi2_gp_model = None

            if not (dataset.kind == "RV" or dataset.kind == "radial_velocity"):
                # TODO remove option 'RV' in PyORBIT version 12
                continue

            if "none" in dataset.models or "None" in dataset.models:
                continue
            if not dataset.models:
                continue

            ndata += dataset.n
            theta_flag = dataset.check_theta_flag(theta_flag)

            for model_name in dataset.models:
                for common_ref in self.models[model_name].common_ref:
                    theta_flag = self.common_models[common_ref].check_theta_flag(
                        theta_flag
                    )

                # TODO: remove try-except starting from PyORBIT version 12 !!
                try:
                    for planet_name in self.models[model_name].multiple_planets:
                        theta_flag = self.common_models[planet_name].check_theta_flag(
                            theta_flag
                        )
                except TypeError:
                    pass

                theta_flag = self.models[model_name].check_theta_flag(
                    dataset_name, theta_flag
                )

                if getattr(self.models[model_name], "internal_likelihood", False):
                    logchi2_gp_model = model_name

            """ Gaussian Process check MUST be the last one or the program will fail
            that's because for the GP to work we need to know the _deterministic_ part of the model
            (i.e. the theoretical values you get when you feed your model with the parameter values) """

            if logchi2_gp_model:
                for common_ref in self.models[logchi2_gp_model].common_ref:
                    theta_flag = self.common_models[common_ref].check_theta_flag(
                        theta_flag
                    )

                theta_flag = self.models[logchi2_gp_model].check_theta_flag(
                    dataset_name, theta_flag
                )

        ndim = np.sum(theta_flag)

        return ndata, ndim, theta_flag

    def rv_log_likelihood(self, theta):

        log_likelihood = 0.00

        if self.dynamical_model is not None:
            """ check if any keyword ahas get the output model from the dynamical tool
            we must do it here because all the planet are involved"""
            dynamical_output = self.dynamical_model.compute(self, theta)

        delayed_lnlk_computation = []
        residuals_analysis = {}

        for dataset_name, dataset in self.dataset_dict.items():
            data_are_rvs = dataset.kind == "RV" or dataset.kind == "radial_velocity"
            # TODO remove option 'RV' in PyORBIT version 12
            logchi2_gp_model = None
            compute_gp_residuals = False

            dataset.model_reset()
            parameter_values = dataset.convert(theta)
            dataset.compute(parameter_values)

            if "none" in dataset.models or "None" in dataset.models:
                continue
            if not dataset.models:
                continue

            skip_loglikelihood = False

            for model_name in dataset.models:
                parameter_values = {}
                for common_ref in self.models[model_name].common_ref:
                    parameter_values.update(
                        self.common_models[common_ref].convert(theta)
                    )

                # TODO: remove try-except starting from PyORBIT version 12 !!
                try:
                    for planet_name in self.models[model_name].multiple_planets:
                        parameter_values.update(
                            self.common_models[planet_name].convert_with_name(
                                theta, planet_name
                            )
                        )
                except TypeError:
                    pass

                parameter_values.update(
                    self.models[model_name].convert(theta, dataset_name)
                )

                """ residuals will be computed following the definition in Dataset class
                    This section has never been tested
                """
                if getattr(self.models[model_name], "residuals_analysis", False):
                    skip_loglikelihood = True
                    if self.models[model_name].gp_before_correlation:
                        compute_gp_residuals = True

                    try:
                        residuals_analysis[model_name]["parameter_values"] = (
                            parameter_values
                        )
                    except:
                        residuals_analysis[model_name] = {
                            "parameter_values": parameter_values
                        }

                    if dataset_name == self.models[model_name].x_dataset:
                        residuals_dataset_label = "x"
                        residuals_analysis[model_name]["x_gp_model"] = None
                        residuals_analysis[model_name]["x_gp_parameters"] = None
                    else:
                        residuals_dataset_label = "y"
                        residuals_analysis[model_name]["y_gp_model"] = None
                        residuals_analysis[model_name]["y_gp_parameters"] = None
                    continue

                if getattr(self.models[model_name], "internal_likelihood", False):
                    logchi2_gp_model = model_name
                    continue

                # if getattr(self.models[model_name], 'model_class', None) is 'common_jitter':
                if getattr(self.models[model_name], "jitter_model", False):
                    dataset.jitter += self.models[model_name].compute(
                        parameter_values, dataset
                    )
                    continue

                if getattr(dataset, "dynamical", False) and getattr(
                    self.models[model_name], "dynamical_model", False
                ):
                    if dataset.kind == "photometry":
                        dataset.buffer_model = self.models[
                            model_name
                        ].compute_dynamical(
                            parameter_values, dataset, dynamical_output[dataset_name]
                        )
                    else:
                        dataset.external_model = dynamical_output[dataset_name]

                if dataset.normalization_model is None and (
                    self.models[model_name].unitary_model
                    or self.models[model_name].normalization_model
                ):
                    dataset.normalization_model = np.ones(dataset.n, dtype=np.double)

                if self.models[model_name].unitary_model:
                    dataset.unitary_model += self.models[model_name].compute(
                        parameter_values, dataset
                    )
                elif self.models[model_name].normalization_model:
                    dataset.normalization_model *= self.models[model_name].compute(
                        parameter_values, dataset
                    )
                else:
                    dataset.additive_model += self.models[model_name].compute(
                        parameter_values, dataset
                    )

            dataset.compute_model()
            dataset.compute_residuals()

            """ Gaussian Process check MUST be the last one or the program will fail
            that's because for the GP to work we need to know the _deterministic_ part of the model
            (i.e. the theoretical values you get when you feed your model with the parameter values) """

            if logchi2_gp_model:
                parameter_values = {}
                for common_ref in self.models[logchi2_gp_model].common_ref:
                    parameter_values.update(
                        self.common_models[common_ref].convert(theta)
                    )

                parameter_values.update(
                    self.models[logchi2_gp_model].convert(theta, dataset_name)
                )

                """ GP Log-likelihood is not computed now because a single matrix must be
                    computed with the joint dataset"""
                if hasattr(self.models[logchi2_gp_model], "delayed_lnlk_computation"):
                    self.models[logchi2_gp_model].add_internal_dataset(
                        parameter_values, dataset
                    )
                    if logchi2_gp_model not in delayed_lnlk_computation:
                        delayed_lnlk_computation.append(logchi2_gp_model)

                elif skip_loglikelihood:
                    residuals_analysis[model_name][
                        residuals_dataset_label + "_gp_model"
                    ] = logchi2_gp_model
                    residuals_analysis[model_name][
                        residuals_dataset_label + "_gp_parameters"
                    ] = parameter_values
                elif data_are_rvs:
                    log_likelihood += self.models[logchi2_gp_model].lnlk_compute(
                        parameter_values, dataset
                    )

                if compute_gp_residuals:
                    dataset.residuals_for_regression -= self.models[
                        logchi2_gp_model
                    ].sample_predict(parameter_values, dataset)
                    residuals_analysis[model_name][
                        residuals_dataset_label + "_gp_model"
                    ] = None

            elif data_are_rvs and (not skip_loglikelihood):
                log_likelihood += dataset.model_logchi2()

        """ Correlation between residuals of two datasets through orthogonal distance regression
            it must be coomputed after any other model has been removed from the
            independent variable, if not provided as ancillary dataset
        """
        for model_name in residuals_analysis:
            parameter_values = residuals_analysis[model_name]["parameter_values"]
            x_dataset = self.dataset_dict[self.models[model_name].x_dataset]
            y_dataset = self.dataset_dict[self.models[model_name].y_dataset]
            modelout_xx, modelout_yy = self.models[model_name].compute(
                parameter_values, x_dataset, y_dataset
            )

            x_dataset.residuals -= modelout_xx
            if residuals_analysis[model_name]["x_gp_model"]:
                logchi2_gp_model = residuals_analysis[model_name]["x_gp_model"]
                parameter_values = residuals_analysis[model_name]["x_gp_parameters"]
                data_log_likelihood = self.models[logchi2_gp_model].lnlk_compute(
                    parameter_values, x_dataset
                )
            else:
                data_log_likelihood += x_dataset.model_logchi2()

            if x_dataset.kind == "RV" or x_dataset.kind == "radial_velocity":
                # TODO remove option 'RV' in PyORBIT version 12
                log_likelihood += data_log_likelihood

            y_dataset.residuals -= modelout_yy
            if residuals_analysis[model_name]["y_gp_model"]:
                logchi2_gp_model = residuals_analysis[model_name]["y_gp_model"]
                parameter_values = residuals_analysis[model_name]["y_gp_parameters"]
                data_log_likelihood = self.models[logchi2_gp_model].lnlk_compute(
                    parameter_values, y_dataset
                )
            else:
                data_log_likelihood = y_dataset.model_logchi2()

            if x_dataset.kind == "RV" or x_dataset.kind == "radial_velocity":
                # TODO remove option 'RV' in PyORBIT version 12
                log_likelihood += data_log_likelihood

        """ In case there is more than one GP model """
        for logchi2_gp_model in delayed_lnlk_computation:
            log_likelihood += self.models[logchi2_gp_model].lnlk_rvonly_compute()

        return log_likelihood

    def log_priors(self, theta):

        log_priors = 0.00

        if not self.include_priors:
            return log_priors

        for model_name, model in self.common_models.items():
            log_priors += model.return_priors(theta)

        for dataset_name, dataset in self.dataset_dict.items():
            dataset.model_reset()
            parameter_values = dataset.convert(theta)
            dataset.compute(parameter_values)

            log_priors += dataset.return_priors(theta)
            if "none" in dataset.models or "None" in dataset.models:
                continue
            if not dataset.models:
                continue
            for model_name in dataset.models:
                log_priors += self.models[model_name].return_priors(theta, dataset_name)

                if getattr(self.models[model_name], "residuals_analysis", False):
                    if dataset_name == self.models[model_name].x_dataset:
                        pass
                    else:
                        # Subtract the priors previously added, to avoid including the priors twice
                        log_priors -= self.models[model_name].return_priors(
                            theta, dataset_name
                        )

        if np.isnan(log_priors):
            log_priors = -np.inf

        return log_priors

    def recenter_bounds(self, pop_mean, recenter=True):
        # This function recenters the bounds limits for circular parameters
        # Also, it extends the range of a parameter if the output of PyDE is a fixed number

        ind_list = []

        for model_name, model in self.common_models.items():
            ind_list.extend(model.index_recenter_bounds())

        for dataset_name, dataset in self.dataset_dict.items():
            for model in dataset.models:
                ind_list.extend(
                    self.models[model].index_recenter_bounds(dataset.name_ref)
                )

        if not recenter:
            return ind_list

        if ind_list:
            tmp_range = (self.bounds[:, 1] - self.bounds[:, 0]) / 2
            replace_bounds = np.zeros([self.ndim, 2])
            replace_bounds[:, 0] = pop_mean - tmp_range
            replace_bounds[:, 1] = pop_mean + tmp_range
            self.bounds[ind_list, :] = replace_bounds[ind_list, :]

    def fix_population(self, pop_mean, population):

        ind_list = self.recenter_bounds(pop_mean, recenter=False)
        n_pop = np.size(population, axis=0)

        if ind_list:
            for par_ind in ind_list:
                fix_sel = (population[:, par_ind] <= self.bounds[par_ind, 0]) | (
                    population[:, par_ind] >= self.bounds[par_ind, 1]
                )
                population[fix_sel, par_ind] = pop_mean[par_ind]

        for ii in range(0, self.ndim):
            if np.amax(population[:, ii]) - np.amin(population[:, ii]) < 10e-14:
                average_pops = np.average(population[:, ii])
                population[:, ii] = np.random.normal(average_pops, 10e-12, n_pop)

        return population
