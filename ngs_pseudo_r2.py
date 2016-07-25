# TC, 7/25/16

import itertools
import numpy as np
import random
import tqdm
# from scipy.misc import logsumexp

#########################
# Signal Allocation Code
#########################

def is_monotonic_list(L):
    """
    Return True if items in list do not decrease, otherwise return False
    """
    return all([(L[m+1] - L[m]) >= 0 for m in range(len(L) - 1)])

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
            if is_monotonic_list(I)]

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

def seq_log_lhood(D, S, L_I, P, prior=None):
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
        prior = np.ones((len(L_I),1), dtype=float) / len(L_I)
    
    # out = logsumexp([seq_alloc_log_lhood(D,S,I,P) for I in L_I], b=prior)
    out = np.max([seq_alloc_log_lhood(D,S,I,P) + np.log(p) for I, p in zip(L_I, prior)])

    return out

#########################
# Population Likelihood Code
#########################

def pop_log_lhood(A_D, S, P):
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
    
    return np.sum(np.apply_along_axis(seq_log_lhood, 1, A_D, S, L_I, P))

def all_log_lhood_fast(A_D, S, P):
    """
    Return a 1-d array of the MAP log-likelihoods of the sequences in A_D
    
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
    possible_Ds = list(itertools.product([0,1], repeat=len(L_D[0])))
    
    prob_Ds = {D: seq_log_lhood(np.array(D),S,L_I,P) for D in possible_Ds}
    
    return np.apply_along_axis(lambda x: prob_Ds[tuple(x)], 1, A_D)
    # return [np.apply_along_axis(lambda x: (x == D).all(), 1, A_D) for D in possible_Ds]
    # return np.array([ct * p for ct, p in zip(count_Ds, prob_Ds)])

def pop_log_lhood_fast(A_D, S, P):
    """
    Return the MAP log-likelihood of a population.

    This is legacy code to work with the slow bootstrapping code.
    """
    return np.sum(all_log_lhood_fast(A_D, S, P))

#########################
# Bootstrap and Pseudo-R2 Code
#########################

def rel_pseudo_r2_pop(L_D, S1, S2, P1, P2=None):
    """
    Return McFadden's Pseudo-R^2 (1-\frac{logL_S1}{logL_S2} 
    of sequences under signal S1 over signal S2.   
    """
    if not P2:
        P2 = P1
    
    logL_1 = pop_log_lhood_fast(L_D, S1, P1)
    logL_2 = pop_log_lhood_fast(L_D, S2, P2)
    
    return 1 - (logL_1/logL_2)

def bootstr_ci(func, X, alpha=0.05, args=None, n_iter=1000, samp_size=None):
    if args:
        f = lambda data: func(data, *args)
    else:
        f = func
        
    if not samp_size:
        samp_size = len(X)
    
    Fstar_lo, Fstar_hi = np.floor(n_iter * np.array([alpha/2., 1-alpha/2.])).astype(int)
    
    bootstrap_stats = np.zeros((n_iter,1))
    
    for i in tqdm.tqdm(range(n_iter)):
        resampled_idx = np.random.randint(X.shape[0], size=(samp_size,1))
        resampled_X = X[resampled_idx,:].squeeze()
        bootstrap_stats[i] = f(resampled_X)
    
    srt_bootstrap_stats = np.sort(bootstrap_stats.squeeze())
    
    ci = (srt_bootstrap_stats[Fstar_lo], srt_bootstrap_stats[Fstar_hi])
    
    return ci

def bootstr_ci_pseudo_r2(L_D, S1, S2, P1, P2=None, alpha=0.05, n_iter=1000, samp_size=None):
    """
    Does bootstrap for pseudo-r2. Gets significant speedups as we can precompute the 
    log-likelihood of all the sequences first, then just bootstrap over the sum.
    """
    if not P2: # This might be better handled by *args...
        P2 = P1
        
    if not samp_size:
        samp_size = L_D.shape[0]
        
    Fstar_lo, Fstar_hi = np.floor(n_iter * np.array([alpha/2., 1-alpha/2.])).astype(int)
    
    # Pre-calculate log-likelihoods
    logLs_1 = all_log_lhood_fast(L_D, S1, P1)
    logLs_2 = all_log_lhood_fast(L_D, S2, P2)
    
    bootstrap_stats = np.zeros((n_iter,1))
    
    for i in tqdm.tqdm(range(n_iter)):
        resampled_idx = np.random.randint(L_D.shape[0], size=(samp_size,1))
        bootstrap_stats[i] = 1 - (np.sum(logLs_1[resampled_idx].squeeze()) /
                                  np.sum(logLs_2[resampled_idx].squeeze()))
    
    srt_bootstrap_stats = np.sort(bootstrap_stats.squeeze())
    
    ci = (srt_bootstrap_stats[Fstar_lo], srt_bootstrap_stats[Fstar_hi])
    return ci