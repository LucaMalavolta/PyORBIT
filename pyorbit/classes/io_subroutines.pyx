import cPickle as pickle
import numpy as np
import copy


def pyde_create_dummy_file(mc,  prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    file_dummy = open(mc.pyde_dir_output + add_prefix + "dummy_file", "wb")
    file_dummy.close()


def pyde_save_to_pickle(mc, population, starting_point, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(mc, open(mc.pyde_dir_output + add_prefix + "model_container.p", "wb"), -1)
    pickle.dump(population, open(mc.pyde_dir_output + add_prefix + "population.p", "wb"), -1)
    pickle.dump(starting_point, open(mc.pyde_dir_output + add_prefix + "starting_point.p", "wb"), -1)


def pyde_load_from_cpickle(pyde_dir_output, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')

    mc = pickle.load(open(pyde_dir_output + add_prefix + "model_container.p", "rb"))
    population = pickle.load(open(pyde_dir_output + add_prefix + "population.p", "rb"))
    starting_point = pickle.load(open(pyde_dir_output + add_prefix + "starting_point.p", "rb"))

    return mc, population, starting_point


def emcee_create_dummy_file(mc, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    file_dummy = open(mc.emcee_dir_output + add_prefix + "dummy_file", "wb")
    file_dummy.close()


def emcee_save_to_cpickle(mc, starting_point, population, prob, state, sampler, samples=None, prefix=None):

    if samples:
        mc.emcee_parameters['nsteps'] = samples
    add_prefix = (prefix + '_' if prefix else '')

    pickle.dump(mc, open(mc.emcee_dir_output + add_prefix + "model_container.p", "wb"), -1)
    pickle.dump(starting_point, open(mc.emcee_dir_output + add_prefix + "starting_point.p", "wb"), -1)
    pickle.dump(population, open(mc.emcee_dir_output + add_prefix + "starting_population.p", "wb"), -1)
    pickle.dump(prob, open(mc.emcee_dir_output + add_prefix + "prob.p", "wb"), -1)
    pickle.dump(state, open(mc.emcee_dir_output + add_prefix + "state.p", "wb"), -1)
    pickle.dump(sampler.chain, open(mc.emcee_dir_output + add_prefix + "sampler_chain.p", "wb"), -1)
    pickle.dump(sampler.lnprobability, open(mc.emcee_dir_output + add_prefix + "sampler_lnprobability.p", "wb"), -1)
    pickle.dump(sampler.acceptance_fraction, open(mc.emcee_dir_output + add_prefix + "sampler_acceptance_fraction.p", "wb"), -1)


def emcee_load_from_cpickle(emcee_dir_output, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')

    mc = pickle.load(open(emcee_dir_output + add_prefix + "model_container.p", "rb"))
    starting_point = pickle.load(open(emcee_dir_output + add_prefix + "starting_point.p", "rb"))
    population = pickle.load(open(emcee_dir_output + add_prefix + "starting_population.p", "rb"))
    prob = pickle.load(open(mc.emcee_dir_output + add_prefix + "prob.p", "rb"))
    state = pickle.load(open(mc.emcee_dir_output + add_prefix + "state.p", "rb"))
    sampler_chain = pickle.load(open(emcee_dir_output + add_prefix + "sampler_chain.p", "rb"))
    sampler_lnprobability = pickle.load(open(emcee_dir_output + add_prefix + "sampler_lnprobability.p", "rb"))
    sampler_acceptance_fraction = pickle.load(open(emcee_dir_output + add_prefix + "sampler_acceptance_fraction.p", "rb"))

    return mc, starting_point, population, prob, state, \
            sampler_chain, sampler_lnprobability, sampler_acceptance_fraction


def nested_sampling_create_dummy_file(mc, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    file_dummy = open(mc.output_directory + add_prefix + "dummy_file", "wb")
    file_dummy.close()


def nested_sampling_save_to_cpickle(mc, prefix=None):

    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(mc, open(mc.output_directory + add_prefix + "model_container.p", "wb"), -1)


def nested_sampling_load_from_cpickle(output_directory, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')
    mc = pickle.load(open(output_directory + add_prefix + "model_container.p", "rb"))
    return mc


def emcee_flatchain(chain, nburnin, nthin):
    """flattening of the emcee chains with removal of burn-in"""
    _, d, _ = np.shape(chain)
    nburn = nburnin / nthin
    if nburn >= d*0.9:
        nburn = d/4

    s = chain[:, nburn:, :].shape
    return chain[:, nburn:, :].reshape(s[0] * s[1], s[2])


def emcee_flatlnprob(lnprob, nburnin, nthin, emcee_version):
    if emcee_version == '3':
        """flattening of the emcee chains with removal of burn-in"""
        d, _ = np.shape(lnprob)
        nburn = nburnin / nthin
        if nburn >= d * 0.9:
            nburn = d / 4

        s = lnprob[nburn:, :].shape
        return lnprob[nburn:, :].reshape(s[0] * s[1])
    else:
        """flattening of the emcee chains with removal of burn-in"""
        _, d = np.shape(lnprob)
        nburn = nburnin / nthin
        if nburn >= d * 0.9:
            nburn = d / 4

        s = lnprob[:, nburn:].shape
        return lnprob[:, nburn:].reshape(s[0] * s[1])


def GelmanRubin(chains_T):
    # Courtesy of Luca "Sbuffo" Borsato
    n, M = np.shape(chains_T)

    theta_m = [np.mean(chains_T[:,i_m]) for i_m in range(0, M)]
    theta = np.mean(theta_m)

    d_theta2 = (theta_m - theta)**2
    B_n = np.sum(d_theta2) / (M-1)

    arg_W = [np.sum((chains_T[:,i_m] - theta_m[i_m])**2) / (n-1) for i_m in range(0, M)]
    W = np.mean(arg_W)

    n_frac = (n-1)/n
    var_plus = n_frac*W + B_n
    Var = var_plus + (B_n/M)

    Rc = np.sqrt(Var / W)
    return Rc

def GelmanRubin_v2(sampler_chain):
    """
    :param chain_T:
    :return:
    """

    """
    from http://joergdietrich.github.io/emcee-convergence.html
    """
    ssq = np.var(sampler_chain, axis=1, ddof=1)
    W = np.mean(ssq, axis=0)
    theta_b = np.mean(sampler_chain, axis=1)
    theta_bb = np.mean(theta_b, axis=0)
    m = sampler_chain.shape[0] * 1.0
    n = sampler_chain.shape[1] * 1.0
    B = n / (m - 1) * np.sum((theta_bb - theta_b)**2, axis=0)
    var_theta = (n - 1) / n * W + 1 / n * B
    Rhat = np.sqrt(var_theta / W)
    return Rhat

"""
def model_container_plot(mc):

    This subroutine makes a deepcopy of the model_container object. Then it substitutes the original datasets with
    densely sampled datasets for plotting purpose

    :param mc:
    :return mc_deepcopy:

    mc_deepcopy = copy.deepcopy(mc)
    mc_deepcopy.deepcopy_for_plot = True
    for dataset_name, dataset in mc_deepcopy.dataset_dict.iteritems():
        if dataset.kind == 'Tcent':
            continue

        x_start = np.min(dataset.x)
        x_end = np.max(dataset.x)
        x_range = x_end-x_start
        x_start -= x_range/20.
        x_end += x_range/20.

        bjd_plot = np.arange(x_start, x_end, 0.1)
        data_input = np.zeros([np.size(bjd_plot, axis=0), 6], dtype=np.double) - 1.
        data_input[:, 0] = bjd_plot[:]

        dataset.define_dataset_base(data_input, update=True)

    mc_deepcopy.model_setup()
    return mc_deepcopy

"""
