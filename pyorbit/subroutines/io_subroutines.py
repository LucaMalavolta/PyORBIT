try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
import numpy as np

try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError


def funcname(self,  parameter_list):
    pass


def pyde_write_dummy_file(mc, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    file_dummy = open(mc.pyde_dir_output + add_prefix + "dummy_file", "wb")
    file_dummy.close()


def pyde_save_to_pickle(mc, population, starting_point, theta_dict, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(mc,
        open(mc.pyde_dir_output + add_prefix + "model_container.p", "wb"))
    pickle.dump(population,
        open(mc.pyde_dir_output + add_prefix + "population.p", "wb"))
    pickle.dump(starting_point,
        open(mc.pyde_dir_output + add_prefix + "starting_point.p", "wb"))
    pickle.dump(theta_dict,
        open(mc.pyde_dir_output + add_prefix + "theta_dictionary.p", "wb"))


def pyde_load_from_cpickle(pyde_dir_output, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')

    mc = pickle.load(
        open(pyde_dir_output + add_prefix + "model_container.p", "rb"))
    population = pickle.load(
        open(pyde_dir_output + add_prefix + "population.p", "rb"))
    starting_point = pickle.load(
        open(pyde_dir_output + add_prefix + "starting_point.p", "rb"))
    theta_dict = pickle.load(
        open(pyde_dir_output + add_prefix + "theta_dictionary.p", "rb"))

    return mc, population, starting_point, theta_dict


def emcee_write_dummy_file(mc, prefix=''):
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
    #pickle.dump(population, open(mc.emcee_dir_output + add_prefix + "starting_population.p", "wb"))
    pickle.dump(population, open(mc.emcee_dir_output + add_prefix + "population.p", "wb"))
    pickle.dump(prob, open(mc.emcee_dir_output + add_prefix + "prob.p", "wb"))
    pickle.dump(state, open(mc.emcee_dir_output + add_prefix + "state.p", "wb"))
    pickle.dump(sampler, open(mc.emcee_dir_output + add_prefix + "sampler.p", "wb"))
    pickle.dump(sampler.chain, open(mc.emcee_dir_output + add_prefix + "sampler_chain.p", "wb"))
    pickle.dump(sampler.lnprobability, open(mc.emcee_dir_output + add_prefix + "sampler_lnprobability.p", "wb"))
    pickle.dump(sampler.acceptance_fraction,
                open(mc.emcee_dir_output + add_prefix + "sampler_acceptance_fraction.p", "wb"))


def zeus_write_dummy_file(mc, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    file_dummy = open(mc.zeus_dir_output + add_prefix + "dummy_file", "wb")
    file_dummy.close()


def zeus_save_to_cpickle(mc, starting_point, population, prob, state, sampler, theta_dict, samples=None, prefix=None):
    if samples:
        mc.zeus_parameters['nsteps'] = samples
    add_prefix = (prefix + '_' if prefix else '')

    pickle.dump(theta_dict, open(mc.zeus_dir_output + add_prefix + "theta_dict.p", "wb"))
    pickle.dump(mc, open(mc.zeus_dir_output + add_prefix + "model_container.p", "wb"))
    pickle.dump(starting_point, open(mc.zeus_dir_output + add_prefix + "starting_point.p", "wb"))
    #pickle.dump(population, open(mc.zeus_dir_output + add_prefix + "starting_population.p", "wb"))
    pickle.dump(population, open(mc.zeus_dir_output + add_prefix + "population.p", "wb"))
    pickle.dump(prob, open(mc.zeus_dir_output + add_prefix + "prob.p", "wb"))
    pickle.dump(state, open(mc.zeus_dir_output + add_prefix + "state.p", "wb"))
    pickle.dump(sampler, open(mc.zeus_dir_output + add_prefix + "sampler.p", "wb"))
    pickle.dump(sampler.chain, open(mc.zeus_dir_output + add_prefix + "sampler_chain.p", "wb"))
    pickle.dump(sampler.lnprobability, open(mc.zeus_dir_output + add_prefix + "sampler_lnprobability.p", "wb"))
    pickle.dump(sampler.acceptance_fraction,
                open(mc.zeus_dir_output + add_prefix + "sampler_acceptance_fraction.p", "wb"))



def affine_load_from_cpickle(dir_output, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    # For backward compatibility
    try:
        theta_dict = pickle.load(open(dir_output + add_prefix + "theta_dict.p", "rb"))
    except FileNotFoundError:
        theta_dict = None

    mc = pickle.load(open(dir_output + add_prefix + "model_container.p", "rb"))
    starting_point = pickle.load(open(dir_output + add_prefix + "starting_point.p", "rb"))

    # For backward compatibility after changing a confusing name
    try:
        population = pickle.load(open(dir_output + add_prefix + "population.p", "rb"))
    except FileNotFoundError:
        population = pickle.load(open(dir_output + add_prefix + "starting_population.p", "rb"))

    prob = pickle.load(open(dir_output + add_prefix + "prob.p", "rb"))
    sampler_chain = pickle.load(open(dir_output + add_prefix + "sampler_chain.p", "rb"))
    sampler_lnprobability = pickle.load(open(dir_output + add_prefix + "sampler_lnprobability.p", "rb"))
    sampler_acceptance_fraction = pickle.load(
        open(dir_output + add_prefix + "sampler_acceptance_fraction.p", "rb"))

    return mc, starting_point, population, prob, \
        sampler_chain, sampler_lnprobability, sampler_acceptance_fraction, theta_dict

def affine_simpler_load_from_cpickle(dir_output, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    try:
        state = pickle.load(open(dir_output + add_prefix + "state.p", "rb"))
        sampler = pickle.load(open(dir_output + add_prefix + "sampler.p", "rb"))
    except (FileNotFoundError, ModuleNotFoundError):
        state = None
        sampler = None

    return state, sampler

def emcee_load_from_cpickle(emcee_dir_output, prefix=''):
    return affine_load_from_cpickle(emcee_dir_output, prefix)

def zeus_simpler_load_from_cpickle(zeus_dir_output, prefix=''):
    return affine_simpler_load_from_cpickle(zeus_dir_output, prefix)

def zeus_load_from_cpickle(zeus_dir_output, prefix=''):
    return affine_load_from_cpickle(zeus_dir_output, prefix)

def emcee_simpler_load_from_cpickle(emcee_dir_output, prefix=''):
    return affine_simpler_load_from_cpickle(emcee_dir_output, prefix)


def starting_point_save_to_cpickle(dir_output, starting_point, bounds, theta_dict, prefix=None):
    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(theta_dict, open(dir_output + add_prefix + "theta_dictionary.p", "wb"))
    pickle.dump(starting_point,
                open(dir_output + add_prefix + "starting_point.p", "wb"))
    pickle.dump(bounds, open(dir_output + add_prefix + "boundaries.p", "wb"))


def starting_point_load_from_cpickle(dir_output, prefix=None):
    add_prefix = (prefix + '_' if prefix else '')
    theta_dict = pickle.load(open(dir_output + add_prefix + "theta_dictionary.p", "rb"))
    starting_point = pickle.load(
        open(dir_output + add_prefix + "starting_point.p", "rb"))
    bounds = pickle.load(open(dir_output + add_prefix + "boundaries.p", "rb"))
    return starting_point, bounds, theta_dict


def nested_sampling_write_dummy_file(mc, prefix=''):
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


def dynesty_results_save_to_cpickle(output_directory, results, prefix=None):
    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(results, open(output_directory + add_prefix + "results.p", "wb"))


def dynesty_results_load_from_cpickle(output_directory, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    results = pickle.load(open(output_directory + add_prefix + "results.p", "rb"))
    return results

def dynesty_results_maxevidence_save_to_cpickle(output_directory, results, prefix=None):
    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(results, open(output_directory + add_prefix + "results_maxevidence.p", "wb"))


def dynesty_results_maxevidence_load_from_cpickle(output_directory, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    results = pickle.load(open(output_directory + add_prefix + "results_maxevidence.p", "rb"))
    return results

def nestle_results_save_to_cpickle(output_directory, results, prefix=None):
    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(results, open(output_directory + add_prefix + "results.p", "wb"))


def nestle_results_load_from_cpickle(output_directory, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    results = pickle.load(open(output_directory + add_prefix + "results.p", "rb"))
    return results

def nestle_results_maxevidence_save_to_cpickle(output_directory, results, prefix=None):
    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(results, open(output_directory + add_prefix + "results_maxevidence.p", "wb"))


def nestle_results_maxevidence_load_from_cpickle(output_directory, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    results = pickle.load(open(output_directory + add_prefix + "results_maxevidence.p", "rb"))
    return results


def ultranest_sampler_save_to_cpickle(output_directory, sampler, prefix=None):
    add_prefix = (prefix + '_' if prefix else '')
    pickle.dump(sampler, open(output_directory + add_prefix + "sampler.p", "wb"))


def ultranest_sampler_load_from_cpickle(output_directory, prefix=''):
    add_prefix = (prefix + '_' if prefix else '')
    sampler = pickle.load(open(output_directory + add_prefix + "sampler.p", "rb"))
    return sampler

def affine_burnin_check(chain, nburnin, nthin, nwalkers=False):
    nburn = int(nburnin / nthin)
    modified = False

    if not nwalkers:
        _, d, _ = np.shape(chain)
    else:
        v1, v2 = np.shape(chain)
        if v1 == nwalkers:
            d = v2
        else:
            d = v1

    if nburn >= d * 0.9:
        nburn = int(d / 4)
        modified = True

    return nburn, modified


def affine_flatchain(chain, nburnin, nthin):
    """flattening of the emcee chains with removal of burn-in"""
    nburn, _ = affine_burnin_check(chain, nburnin, nthin)
    s = chain[:, nburn:, :].shape
    return chain[:, nburn:, :].reshape(s[0] * s[1], s[2])


def affine_flatlnprob(lnprob, nburnin, nthin, population, nwalkers):

    nburn, _  = affine_burnin_check(lnprob, nburnin, nthin, nwalkers)

    v1, v2 = np.shape(lnprob)
    if v1 == nwalkers:
        s = lnprob[:, nburn:].shape
        return lnprob[:, nburn:].reshape(s[0] * s[1]), lnprob.T
    else:
        s = lnprob[nburn:, :].shape
        return lnprob[nburn:, :].reshape(s[0] * s[1]), lnprob


def emcee_burnin_check(chain, nburnin, nthin, nwalkers=False):
    return affine_burnin_check(chain, nburnin, nthin, nwalkers)

def zeus_burnin_check(chain, nburnin, nthin, nwalkers=False):
    return affine_burnin_check(chain, nburnin, nthin, nwalkers)

def emcee_flatchain(chain, nburnin, nthin):
    return affine_flatchain(chain, nburnin, nthin)

def emcee_flatlnprob(lnprob, nburnin, nthin, population, nwalkers):
    return affine_flatlnprob(lnprob, nburnin, nthin, population, nwalkers)

def zeus_flatchain(chain, nburnin, nthin):
    return affine_flatchain(chain, nburnin, nthin)

def zeus_flatlnprob(lnprob, nburnin, nthin, population, nwalkers):
    return affine_flatlnprob(lnprob, nburnin, nthin, population, nwalkers)


def GelmanRubin(chains_T):
    # Courtesy of Luca "Sbuffo" Borsato
    n, M = np.shape(chains_T)

    theta_m = [np.mean(chains_T[:, i_m]) for i_m in range(0, M)]
    theta = np.mean(theta_m)

    d_theta2 = (theta_m - theta) ** 2
    B_n = np.sum(d_theta2) / (M - 1)

    arg_W = [np.sum((chains_T[:, i_m] - theta_m[i_m]) ** 2) / (n - 1) for i_m in range(0, M)]
    W = np.mean(arg_W)

    n_frac = (n - 1) / n
    var_plus = n_frac * W + B_n
    Var = var_plus + (B_n / M)

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
    B = n / (m - 1) * np.sum((theta_bb - theta_b) ** 2, axis=0)
    var_theta = (n - 1) / n * W + 1 / n * B
    Rhat = np.sqrt(var_theta / W)
    return Rhat
