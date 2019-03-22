try:
    import cPickle as pickle
except:
    import pickle
import numpy as np


def pyde_create_dummy_file(mc, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    file_dummy = open(mc.pyde_dir_output + add_prefix + "dummy_file", "wb")
    file_dummy.close()


def pyde_save_to_pickle(mc, population, starting_point, theta_dict, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(mc, open(mc.pyde_dir_output + add_prefix + "model_container.p", "wb"))
    pickle.dump(population, open(mc.pyde_dir_output + add_prefix + "population.p", "wb"))
    pickle.dump(starting_point, open(mc.pyde_dir_output + add_prefix + "starting_point.p", "wb"))
    pickle.dump(theta_dict, open(mc.pyde_dir_output + add_prefix + "theta_dictionary.p", "wb"))


def pyde_load_from_cpickle(pyde_dir_output, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')

    mc = pickle.load(open(pyde_dir_output + add_prefix + "model_container.p", "rb"))
    population = pickle.load(open(pyde_dir_output + add_prefix + "population.p", "rb"))
    starting_point = pickle.load(open(pyde_dir_output + add_prefix + "starting_point.p", "rb"))
    theta_dict = pickle.load(open(pyde_dir_output + add_prefix + "theta_dictionary.p", "rb"))

    return mc, population, starting_point, theta_dict


def emcee_create_dummy_file(mc, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    file_dummy = open(mc.emcee_dir_output + add_prefix + "dummy_file", "wb")
    file_dummy.close()


def emcee_save_to_cpickle(mc, starting_point, population, prob, state, sampler, theta_dict, samples=None, prefix=None):

    if samples:
        mc.emcee_parameters['nsteps'] = samples
    add_prefix = (prefix + '_' if prefix else '')

    pickle.dump(theta_dict, open(mc.emcee_dir_output + add_prefix + "theta_dict.p", "wb"))
    pickle.dump(mc, open(mc.emcee_dir_output + add_prefix + "model_container.p", "wb"))
    pickle.dump(starting_point, open(mc.emcee_dir_output + add_prefix + "starting_point.p", "wb"))
    pickle.dump(population, open(mc.emcee_dir_output + add_prefix + "starting_population.p", "wb"))
    pickle.dump(prob, open(mc.emcee_dir_output + add_prefix + "prob.p", "wb"))
    pickle.dump(state, open(mc.emcee_dir_output + add_prefix + "state.p", "wb"))
    pickle.dump(sampler.chain, open(mc.emcee_dir_output + add_prefix + "sampler_chain.p", "wb"))
    pickle.dump(sampler.lnprobability, open(mc.emcee_dir_output + add_prefix + "sampler_lnprobability.p", "wb"))
    pickle.dump(sampler.acceptance_fraction, open(mc.emcee_dir_output + add_prefix + "sampler_acceptance_fraction.p", "wb"))


def emcee_load_from_cpickle(emcee_dir_output, prefix=''):

    add_prefix = (prefix + '_' if prefix else '')

    # For backward compatibility
    try:
        theta_dict = pickle.load(open(emcee_dir_output + add_prefix + "theta_dict.p", "rb"))
    except:
        theta_dict = None

    mc = pickle.load(open(emcee_dir_output + add_prefix + "model_container.p", "rb"))
    starting_point = pickle.load(open(emcee_dir_output + add_prefix + "starting_point.p", "rb"))
    population = pickle.load(open(emcee_dir_output + add_prefix + "starting_population.p", "rb"))
    prob = pickle.load(open(mc.emcee_dir_output + add_prefix + "prob.p", "rb"))
    state = pickle.load(open(mc.emcee_dir_output + add_prefix + "state.p", "rb"))
    sampler_chain = pickle.load(open(emcee_dir_output + add_prefix + "sampler_chain.p", "rb"))
    sampler_lnprobability = pickle.load(open(emcee_dir_output + add_prefix + "sampler_lnprobability.p", "rb"))
    sampler_acceptance_fraction = pickle.load(open(emcee_dir_output + add_prefix + "sampler_acceptance_fraction.p", "rb"))

    return mc, starting_point, population, prob, state, \
            sampler_chain, sampler_lnprobability, sampler_acceptance_fraction, theta_dict


def starting_point_save_to_cpickle(dir_output, starting_point, theta_dict, prefix=None):

    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(theta_dict, open(dir_output + add_prefix + "theta_dictionary.p", "wb"))
    pickle.dump(starting_point, open(dir_output + add_prefix + "starting_point.p", "wb"))


def starting_point_load_from_cpickle(dir_output, prefix=None):

    add_prefix = (prefix + '_' if prefix else '')
    theta_dict = pickle.load(open(dir_output + add_prefix + "theta_dictionary.p", "rb"))
    starting_point = pickle.load(open(dir_output + add_prefix + "starting_point.p", "rb"))
    return starting_point, theta_dict


def nested_sampling_create_dummy_file(mc, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    file_dummy = open(mc.output_directory + add_prefix + "dummy_file", "wb")
    file_dummy.close()


def nested_sampling_save_to_cpickle(mc, prefix=None):

    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(mc, open(mc.output_directory + add_prefix + "model_container.p", "wb"))


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
