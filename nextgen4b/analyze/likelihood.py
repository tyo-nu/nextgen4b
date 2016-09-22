# TC, 7/25/16

import itertools
import numpy as np
import random
import tqdm
from scipy.misc import logsumexp

#########################
# Bootstrap and Pseudo-R2 Code
#########################

def load_words_to_array(fname, rare_base='A',
                        discard_base='-', del_as_misinc=False):
    """
    Return a KxM {0,1} array representing the M-length words contained in
    one of the word.txt files you've been using
    """
    with open(fname) as i_f:
        txt_data = [line.strip() for line in i_f.readlines()
                    if discard_base not in line]
    
    if del_as_misinc:
        rare_base_list = [rare_base, 'd']
    else:
        rare_base_list = [rare_base]
    return np.array([np.where([c in rare_base_list for c in s], 0, 1)
                     for s in txt_data])

def load_controls_to_P(A_D_0, A_D_1):
    p1 = np.sum(A_D_0, axis=0) / float(A_D_0.shape[0])
    p2 = np.sum(A_D_1, axis=0) / float(A_D_1.shape[0])

    return np.array([p1, p2]).T

#########################
# Signal Allocation Code
#########################

def is_monotonic(vec):
    """
    Return True if items in list do not decrease, otherwise return False
    """
    return all([(vec[m+1] - vec[m]) >= 0 for m in range(len(vec) - 1)])

def get_monotonic_allocations(length_d, length_s):
    """
    Return a list of possible ways the strand D could have be written during some signal S.

    Inputs:
        length_d - The length of the DNA analysis strand
        length_s - The number of signal "phases" 
    Output:
        I_L - A list of lists, I. Individual I's are of length D, with each element being
              an index of an element of S, indicating that the nucleotide was written during
              that portion of the signal.
    """

    return [list(I) for I in itertools.product(range(length_s), repeat=length_d)
            if is_monotonic(I)]

#########################
# Single-Sequence Likelihood Code
#########################

def seq_alloc_log_lhood(D, S, I, P):
    """
    Get the likelihood of strand D under signal S, given positional probabilities P. I
    is an "allocation" that indexes D to S.

    Input:
        D - A 1-d M-length array composed of {0,1}, representing errors on a DNA strand
        S - A 1-d N-length array composed of {0,1}, representing the signal delivered at each phase of the signal
        I - An M-length list, providing indexing of each element of D to an element of S
        P - A numpy Nx2 array, each row composed of error rates [p_{m,0}, p_{m,1}]

    Output:
        ll - The log-likelihood of the strand given S, I, and P
    """
    # Get the vector of error rates
    # First get the right S_n index for the d_m's
    # Then use where to get the right probabilities from P
    ers = np.where(S[I], P[:, 1], P[:, 0])

    lhood_components = D*ers + (1-D)*(1-ers)
    return np.sum(np.log(lhood_components))

def seq_log_lhood(D, S, L_I, P, prior=None, method='map'):
    """
    Return the MAP log-likelihood of a sequence over all possible allocations of that
    sequence over the signal.
    
    Input:
        D - A numpy Mx1 array composed of {0,1}, representing errors on a DNA strand
        S - A numpy Nx1 array composed of {0,1}, representing the signal delivered at each phase of the signal
        L_I - A list of M-length lists, providing indexing of each element of D to an element of S
        P - A numpy Nx2 array, each row composed of error rates [p_{m,0}, p_{m,1}]
        prior - A numpy array the length of L_I. The likelihood of each allocation, defaults to a uniform prior

    Output:
        ll - The maximum log-likelihood of the strand given S, I, and P
    """
    if not prior: # i.e. it's uniform
        prior = np.ones(len(L_I), dtype=float) / len(L_I)

    if method is 'map':
        out = np.max([seq_alloc_log_lhood(D,S,I,P) + np.log(p)
                      for I, p in zip(L_I, prior)])
    elif method is 'sum':
        out = logsumexp([seq_alloc_log_lhood(D,S,I,P) for I in L_I], b=prior)
    else:
        raise ValueError('Unrecognized method `%s`' % method)

    return out

#########################
# Population Likelihood Code
#########################

def pop_log_lhood(A_D, S, P, **kwargs):
    """
    Return the log-likelihood of a population of sequences.
    
    Input:
        A_D - A numpy KxM array composed of {0,1}, representing errors on a M-length DNA strand for K strands
        S - A numpy Nx1 array composed of {0,1}, representing the signal delivered at each phase of the signal
        P - A numpy Nx2 array, each row composed of error rates [p_{m,0}, p_{m,1}]
        
    Output:
        ll - The log-likelihood of the population given S and P
    """
    
    assert A_D.size
    L_I = get_monotonic_allocations(A_D.shape[1], len(S))
    
    return np.sum(np.apply_along_axis(seq_log_lhood, 1, A_D, S, L_I, P, **kwargs))

def all_log_lhood_fast(A_D, S, P, **kwargs):
    """
    Return a 1-d array of the log-likelihoods of the sequences in A_D
    
    Input:
        A_D - A numpy KxM array composed of {0,1}, representing errors on a M-length DNA strand for K strands
        S - A 1-d N-length array composed of {0,1}, representing the signal delivered at each phase of the signal
        P - A numpy Nx2 array, each row composed of error rates [p_{m,0}, p_{m,1}]
        
    Output:
        ll - The log-likelihood of the population given S and P
    """
    assert A_D.size
    L_I = get_monotonic_allocations(A_D.shape[1], len(S))
    
    # Get all possible D's
    possible_Ds = list(itertools.product([0,1], repeat=len(A_D[0])))
    
    prob_Ds = {D: seq_log_lhood(np.array(D), S, L_I, P, **kwargs)
               for D in possible_Ds}
    
    return np.apply_along_axis(lambda x: prob_Ds[tuple(x)], 1, A_D)
    # return [np.apply_along_axis(lambda x: (x == D).all(), 1, A_D) for D in possible_Ds]
    # return np.array([ct * p for ct, p in zip(count_Ds, prob_Ds)])

def pop_log_lhood_fast(A_D, S, P, **kwargs):
    """
    Return the log-likelihood of a population, using a method in kwargs.
    """
    return np.sum(all_log_lhood_fast(A_D, S, P, **kwargs))

def multi_pop_log_lhood_fast(l_A_D, S, P, **kwargs):
    """
    Return the log-likelihood of a multi-sample population.
    """
    return np.sum([np.sum(all_log_lhood_fast(A_D, S, P, **kwargs))
                   for A_D in l_A_D])

#########################
# Bootstrap and Pseudo-R2 Code
#########################

def rel_pseudo_r2_pop(L_D, S1, S2, P1, P2=None, **kwargs):
    """
    Return McFadden's Pseudo-R^2 (1-\frac{logL_S1}{logL_S2} 
    of sequences under signal S1 over signal S2.   
    """
    if not P2:
        P2 = P1
    
    logL_1 = pop_log_lhood_fast(L_D, S1, P1, **kwargs)
    logL_2 = pop_log_lhood_fast(L_D, S2, P2, **kwargs)
    
    return 1 - (logL_1/logL_2)

def multi_rel_pseudo_r2_pop(L_D_list, S1, S2, P1, P2=None, **kwargs):
    """
    Return McFadden's Pseudo-R^2 (1-\frac{logL_S1}{logL_S2} 
    of sequences under signal S1 over signal S2.   
    """
    if not P2:
        P2 = P1
    
    logL_1 = multi_pop_log_lhood_fast(L_D_list, S1, P1, **kwargs)
    logL_2 = multi_pop_log_lhood_fast(L_D_list, S2, P2, **kwargs)
    
    return 1 - (logL_1/logL_2)

def bootstr_ci(func, X, alpha=0.05, args=None, kwargs=None,
               n_iter=1000, samp_size=None, verbose=False):
    if args and kwargs:
        f = lambda data: func(data, *args, **kwargs)
    elif args:
        f = lambda data: func(data, *args)
    elif kwargs:
        f = lambda data: func(data, **kwargs)
    else:
        f = func

    if not samp_size:
        samp_size = len(X)
    
    Fstar_lo, Fstar_hi = np.floor(n_iter * np.array([alpha/2., 1-alpha/2.])).astype(int)

    f_mu = f(X)
    
    bootstrap_stats = np.zeros((n_iter,) + f_mu.shape)
    
    if verbose:
        for i in tqdm.tqdm(range(n_iter)):
            resampled_idx = np.random.randint(X.shape[0], size=(samp_size,1))
            resampled_X = X[resampled_idx,...].squeeze()
            bootstrap_stats[i,...] = f(resampled_X)
    else:
        for i in range(n_iter):
            resampled_idx = np.random.randint(X.shape[0], size=(samp_size,1))
            resampled_X = X[resampled_idx,...].squeeze()
            bootstrap_stats[i,...] = f(resampled_X)

    srt_bootstrap_stats = np.sort(bootstrap_stats, axis=0)
    
    ci = (srt_bootstrap_stats[Fstar_lo], srt_bootstrap_stats[Fstar_hi])
    
    return ci

def bootstr_ci_pseudo_r2(L_D, S1, S2, P1, P2=None,
                         alpha=0.05, n_iter=1000, samp_size=None,
                         verbose=False, **kwargs):
    """
    Does bootstrap for pseudo-r2. Gets significant speedups as we can precompute the 
    log-likelihood of all the sequences first, then just bootstrap over the sum.
    """
    if not P2: # This might be better handled by *args...
        P2 = P1
    if not samp_size:
        samp_size = L_D.shape[0]

    if verbose:
        iterator = tqdm.tqdm(range(n_iter))
    else:
        iterator = range(n_iter)

    Fstar_lo, Fstar_hi = np.floor(n_iter * np.array([alpha/2., 1-alpha/2.])).astype(int)

    # Pre-calculate log-likelihoods
    logLs_1 = all_log_lhood_fast(L_D, S1, P1, **kwargs)
    logLs_2 = all_log_lhood_fast(L_D, S2, P2, **kwargs)

    bootstrap_stats = np.zeros((n_iter))

    for i in iterator:
        resampled_idx = np.random.randint(L_D.shape[0], size=(samp_size,1))
        bootstrap_stats[i] = 1 - (np.sum(logLs_1[resampled_idx].squeeze()) /
                                  np.sum(logLs_2[resampled_idx].squeeze()))

    srt_bootstrap_stats = np.sort(bootstrap_stats.squeeze())

    ci = (srt_bootstrap_stats[Fstar_lo], srt_bootstrap_stats[Fstar_hi])
    return ci

def multi_bootstr_ci_pseudo_r2(L_D_list, S1, S2, P1, P2=None,
                               alpha=0.05, n_iter=1000,
                               verbose=False, **kwargs):
    """
    Does bootstrap for pseudo-r2. Gets significant speedups as we can
    precompute the log-likelihood of all the sequences first, then just
    bootstrap over the sum.

    Inputs:
    - L_D_list: list of arrays denoting errors for each strand in each sample
    - S1: array denoting the first signal
    - S2: array denoting the second signal
    - P1: array denoting the error rates at each site for each condition in S1
            or S2

    Output:
    - ci: a tuple containing the lower and upper pseudo-R2 bound
    """

    if not P2: # This might be better handled by *args...
        P2 = P1

    if verbose:
        iterator = tqdm.tqdm(range(n_iter))
    else:
        iterator = range(n_iter)

    Fstar_lo, Fstar_hi = np.floor(n_iter * np.array([alpha/2., 1-alpha/2.])).astype(int)

    # Pre-calculate log-likelihoods
    logLs_1 = [all_log_lhood_fast(L_D, S1, P1, **kwargs) for L_D in L_D_list]
    logLs_2 = [all_log_lhood_fast(L_D, S2, P2, **kwargs) for L_D in L_D_list]

    bootstrap_stats = np.zeros((n_iter))

    for i in iterator:
        boot_idx = [np.random.randint(L_D.shape[0], size=(L_D.shape[0], 1))
                    for L_D in L_D_list]

        boot_L1 = [np.sum(logL[bix].squeeze())
                   for logL, bix in zip(logLs_1, boot_idx)]
        boot_L2 = [np.sum(logL[bix].squeeze())
                   for logL, bix in zip(logLs_2, boot_idx)]

        bootstrap_stats[i] = 1 - (np.sum(boot_L1) / np.sum(boot_L2))

    srt_bootstrap_stats = np.sort(bootstrap_stats.squeeze())

    ci = (srt_bootstrap_stats[Fstar_lo], srt_bootstrap_stats[Fstar_hi])
    return ci

def model_softmax(L_D, S_list, P1, **kwargs):
    """
    Calculates softmax of models based on model likelihood
    """
    x = np.array([pop_log_lhood_fast(L_D, S, P1, **kwargs) for S in S_List])

    return np.exp(x) / np.sum(np.exp(x), axis=0)