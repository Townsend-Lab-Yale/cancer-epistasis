"""Cancer epistasis -- Jorge A. Alfaro-Murillo

This module estimates fluxes from and to combinations of mutations.

The total number of mutations can be set with the parameter n.

A combination of mutations is represented by a vector of 0's and 1's,
that has a 1 or a 0 in the i-th position if the individual has or does
not have the i-th mutation.

The likelihood is a Multinomial distribution where the parameters are
given by the relationship between the fluxes and the probability of
being at each state at time t=T, where at t=0 everyone is starting
without mutations. The observed variables are the samples, and the
priors for the fluxes are uninformative uniform distributions.

"""

import os
import numpy as np

import pymc3 as pm
import theano.tensor as tt
from theano.compile.ops import as_op

from theory import numbers_positive_lambdas
from theory import build_S_as_array
from theory import obtain_pos_lambdas_indices
from theory import order_pos_lambdas

from scipy.stats import chi2

random_seed = 777
"""Random seed to feed the random generators, to be able to replicate
results."""
np.random.RandomState(np.random.MT19937(np.random.SeedSequence(random_seed)))


T = 1
"""Time at which the probabilities are finally evaluated"""


## We will compute numerically integrals according to this resolution
resolution = 10000
times = np.linspace(0, T, resolution)


## * Methods
## ** Computing methods

def establish_non_definable_lambdas(samples):
    """Establish which lambdas are actually definable. If the samples of the
    combination mutations x and y are both 0, then the flux from x
    to y is not definable.

    :type samples: numpy.ndarray
    :param samples: A one dimensional array of size 2^M containing the
        number of individuals at each state of S.

    :rtype: numpy.ndarray
    :return: One dimensional array of the same size as :const`S` with
        the computed probilities.

    """
    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations

    S = build_S_as_array(M)

    non_definable = [[True
                      if np.sum(S[j] >= S[i]) == M
                      and np.sum(S[j] - S[i]) == 1
                      and i < len(S)-1 # until here positive lambdas
                      and samples[i] == 0
                      and samples[j] == 0
                      else False
                      for i in range(len(S))]
                     for j in range(len(S))]
    return np.array(non_definable)


def compute_Ps_at_T(positive_lambdas):
    """Compute for each state in S the probability of being at that state
    at time t=T if starting without mutations.

    :type positive_lambdas: numpy.ndarray
    :param positive_lambdas: One dimensional array with the lambda's
        that are positive. It should be of size equal to the number of
        True indeces in `const:positive_lambdas_indices`.

    :rtype: numpy.ndarray
    :return: One dimensional array of the same size as :const`S` with
        the computed probilities.

    """

    M = numbers_positive_lambdas.index(len(positive_lambdas))

    S = build_S_as_array(M)

    lambdas = np.zeros([len(S), len(S)])

    positive_lambdas_indices = obtain_pos_lambdas_indices(S)

    lambdas[positive_lambdas_indices] = positive_lambdas

    lambdas[np.eye(len(S), len(S), dtype=bool)] = -np.sum(lambdas, axis=0)

    Ps = np.zeros((len(S), resolution))

    Ps[0] = np.exp(lambdas[0, 0]*times)

    mutations_per_x = np.sum(S, axis=1)

    lambda_xx = lambdas[0, 0]
    lambda_xys = lambdas[mutations_per_x == 1, 0]
    lambda_yys = lambdas[mutations_per_x == 1, mutations_per_x == 1]
    Ps[mutations_per_x == 1] = (
        lambda_xys[:, np.newaxis] / (lambda_xx - lambda_yys[:, np.newaxis])
        * (np.exp(lambda_xx*times) - np.exp(lambda_yys[:, np.newaxis]*times)))

    for m in range(2, M):
        for k, y in enumerate(S):
            if np.sum(y) == m:
                y_minus_indices = np.array(
                    [True
                     if np.sum(y >= S[i]) == M
                     and np.sum(y - S[i]) == 1
                     else False
                     for i in range(len(S))])
                lambda_yy = lambdas[k, k]
                lambda_y_minus_eis_y = lambdas[k, y_minus_indices]
                ## Notice that Ps[y_minus_indices] were already computed
                to_int = np.exp(-lambda_yy*times) * Ps[y_minus_indices]
                ## trapezoidal rule for integration
                integral = ((np.cumsum(to_int, axis=-1) - to_int/2)
                            * T/resolution)
                Ps[k] = np.sum((lambda_y_minus_eis_y[:, np.newaxis]
                                * np.exp(lambda_yy*times)
                                * integral),
                               axis=0)

    ## This saves time on the last calculation, as we know that the sum
    ## of P_x(t) on x should be 1
    Ps[-1] = 1 - np.sum(Ps[:-1], axis=0)

    return Ps[:, -1]


## as_op is a decorator that transforms functions so that they can be
## used in pymc3 models. It could use the syntatic sugar:
##
## @as_op
## compute_Ps_at_T
##
## and this would have redefined compute_Ps_at_T, but I'll name a new
## function to use the normal compute_Ps_at_T as well
compute_Ps_at_T_tens = as_op(itypes=[tt.dvector],
                             otypes=[tt.dvector])(compute_Ps_at_T)


## ** Main method

def estimate_lambdas(samples, upper_bound_prior=3, draws=10000,
                     burn=500, chains=8, save_name=None, kwargs=None):
    """Run the main simulation to find the estimates of the fluxes.

    If draws is equal to 1 find the maximum a posteriori (MAP)
    estimate, else perform a Metropolis-Hastings Markov chain Monte
    Carlo to evaluate the posterior of the fluxes.

    The fluxes are assumed to have a non-informative [0,1] uniform as
    prior.

    :type samples: numpy.ndarray
    :param samples: One dimensional array with the samples. It should
        be of the same size as S, and have in each entry the number of
        individuals that have the respective mutation combination.

    :type upper_bound_prior: float
    :param upper_bound_prior: Upper bound for the uniform prior used
        for the each lambda.

    :type draws: int
    :param draws: Number of samples to draw.

    :type burn: int
    :param burn: Number of samples to be discarded as the procedure is
        tuning in. Only takes effect if `draws` > 1.

    :type chains: int
    :param chains: Number of chains to run in parallel. Only takes
        effect if `draws` > 1.

    :type kwargs: dict
    :param kwargs: Dictionary of keyword arguments to pass to the pymc3
        find_MAP function if `draws` = 1, or to the sample function if
        `draws` > 1.

    :type save_name: str or NoneType
    :param save_name: Save results under this name. If None (default),
        do not save the results. If `draws` = 1, that is, if the
        results are the MAP estimate, then the dictionary with the
        results is saved as a npz file. If `draws` > 1, then each
        chain of the trace goes inside a directory with this name.

    :rtype: pymc3.backends.base.MultiTrace or dict or tuple
    :return: Estimates of the fluxes as a MultiTrace pymc3 object or,
        if `draws`=1, a dictionary with the MAP (also MLE, because the
        priors are uniform) estimates, or if `draws`=1 and
        'return_raw' is passed to kwargs as True, a tuple with the
        full output of scipy.optimize.minimize.

    """

    M = int(np.log2(len(samples))) # 2^M is the number mutation combinations

    S = build_S_as_array(M)

    positive_lambdas_indices = obtain_pos_lambdas_indices(S)

    number_positive_lambdas = np.sum(positive_lambdas_indices)

    number_samples = np.sum(samples)

    with pm.Model():

        ## We set an uninformative prior for the lambdas:
        positive_lambdas = pm.Uniform(
            name="lambdas",
            lower=0,
            upper=upper_bound_prior,
            shape=number_positive_lambdas)

        Ps = compute_Ps_at_T_tens(positive_lambdas)

        likelihood = pm.Multinomial(name="samples",
                                    p=Ps,
                                    n=number_samples,
                                    observed=samples)
        if kwargs is None:
           kwargs = {}

        if draws == 1:
            results = pm.find_MAP(**kwargs)

        else:
            results = pm.sample(int(draws/chains),
                                cores=chains,
                                tune=burn,
                                step=pm.Metropolis(),
                                random_seed=random_seed,
                                **kwargs)

    if save_name is not None:
        if draws == 1:
            np.savez(save_name, **results)
        else:
            pm.save_trace(results, save_name, overwrite=True)

    return results


## ** Other methods to handle the results

def convert_lambdas_to_dict(results):
    """Convert a fluxes array to a dictionary.

    The indexes of the dictionary are the subscripts of the
    lambdas. So if lambda[x, y] represents the flux from x to y, it
    can be obtain from the result of this function easily by using
    output[x, y], with x and y written as tuples.

    :type results: pymc3.backends.base.MultiTrace or dict or list
    :param results: Estimates of the fluxes as a MultiTrace pymc3
        object or a dictionary with the ordered MAP estimates (both as
        returned by :func:`estimate_lambdas`), or as a list of
        confidence intervals for the fluxes as returned by
        asymp_CI_lambdas.

    :rtype: dict
    :return: Dictionary with the fluxes, indexed by a pair of tuples
        representing the mutation combination where the flux is coming
        from and going to.

    """
    if isinstance(results, dict):
        # if results is a dictionary, it must be the return from
        # estimate_lambdas with draws=1, that is, the MAP (and MLE)
        results = results['lambdas']
    elif isinstance(results, pm.backends.base.MultiTrace):
        # if results is a pymc3.backends.base.MultiTrace, it must be
        # the return from estimate_lambdas with draws>1, that is, the
        # posterior estimates
        results = results.get_values('lambdas').T
    elif not isinstance(results, list):
        # otherwise results should be the confidence intervales and
        # thus should be a list because that is the return of
        # asymp_CI_lambdas
        raise Exception("`results` type not compatible")

    M = numbers_positive_lambdas.index(len(results))

    S = build_S_as_array(M)

    subscripts = order_pos_lambdas(S)

    results_as_dict = {subscript_pair:value
                       for subscript_pair, value in zip(subscripts,
                                                        results)}

    return results_as_dict


def compute_gammas(lambdas, mus):
    """Compute the selection coefficients.

    :type lambdas: dict
    :param lambdas: Dictionary with the fluxes, indexed by a pair of
        tuples representing the mutation combination where the flux is
        coming from and going to. Each flux can be represented by a
        single estimate or by confidence interval given as a list.

    :type mus: dict
    :param mus: Dictionary with the mutation rates indexed by a tuple
        of 0's and one 1, representing which mutation the respective
        rate represents. It is assumed that the mutation rates do not
        change via epistatic effects

    :rtype: dict
    :return: Dictionary with the selection coefficients, indexed by a
        pair of tuples representing the mutation combination where the
        flux is coming from and going to.

    """
    gammas = {
        xy:([flux[0]/mus[tuple(np.array(xy[1])-np.array(xy[0]))],
             flux[1]/mus[tuple(np.array(xy[1])-np.array(xy[0]))]]
            if isinstance(flux, list) else
            flux/mus[tuple(np.array(xy[1])-np.array(xy[0]))])
        for xy, flux in lambdas.items()}
    return gammas


def compute_log_lh(positive_lambdas, samples):
    """Not really the log likelihood but that plus a constant (that
    is, the log of number that is proportional to the likelihood).

    :type positive_lambdas: numpy.ndarray
    :param positive_lambdas: One dimensional array with the lambda's
        that are positive. It should be of size equal to the number of
        True indeces in `const:positive_lambdas_indices`.

    :type samples: numpy.ndarray
    :param samples: A one dimensional array of size 2^M containing the
        number of individuals at each state of S.

    :rtype:
    :return: Dictionary with the fluxes, indexed by a pair of tuples
        representing the mutation combination where the flux is coming
        from and going to.

    """
    Ps = compute_Ps_at_T(positive_lambdas)
    return np.sum(samples*np.log(Ps))


## ** Methods to estimate confidence intervals

def asymp_CI_lambda(lambda_index, lambdas_mle, samples, ci=0.95, tol=10**(-6)):
    """Compute an asymptomatic confidence interval for the flux indexed by
    `lambda_index`.

    An asymptotic confidence interval for a particular lambda is given by:

    {lambda | log_lh(lambda, rest of lambdas_mle) >= log_lh(lambdas_mle) - Chi^2(1,ci)/2}

    where Chi^2(1,ci) is the ci-quantile of the Chi-squared distribution with 1 dof.

    (proof at https://statproofbook.github.io/P/ci-wilks.html)

    """
    rhs = (compute_log_lh(lambdas_mle, samples) - chi2.ppf(ci, 1)/2)

    upper_est = lambdas_mle[lambda_index]
    to_compare = 2*upper_est

    lambdas = lambdas_mle.copy()
    lambdas[lambda_index] = to_compare

    diff = compute_log_lh(lambdas, samples) - rhs

    while np.abs(diff) > tol:
        if diff > 0:
            ## lambdas in CI
            upper_est = to_compare
            to_compare = 2*to_compare
        else:
            ## lambdas not in CI
            to_compare = (to_compare + upper_est)/2

        lambdas[lambda_index] = to_compare
        diff = compute_log_lh(lambdas, samples) - rhs

    upper_est = to_compare


    lower_est = lambdas_mle[lambda_index]
    if lower_est > tol:
        ## otherwise it's practically 0
        to_compare = lower_est/2

        lambdas[lambda_index] = to_compare

        diff = compute_log_lh(lambdas, samples) - rhs

        while np.abs(diff) > tol and to_compare > tol:
            if diff > 0:
                ## lambdas in CI
                lower_est = to_compare
                to_compare = to_compare/2
            else:
                ## lambdas not in CI
                to_compare = (to_compare + lower_est)/2

            lambdas[lambda_index] = to_compare
            diff = compute_log_lh(lambdas, samples) - rhs

        lower_est = to_compare

    return [lower_est, upper_est]


def asymp_CI_lambdas(lambdas_mle, samples, ci=0.95, tol=10**(-6),
                     print_progress=False):

    CIs = []

    for lambda_index in range(len(lambdas_mle)):
        if print_progress:
            print(f"computing CI for flux {lambda_index}/{len(lambdas_mle)}")
        CIs.append(asymp_CI_lambda(lambda_index,
                                   lambdas_mle,
                                   samples,
                                   ci=ci,
                                   tol=tol))

    return CIs
