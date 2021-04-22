#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 22 23:19:41 2021

@author: ptruong
"""

import numpy as np
import numpy.random as npr
import pandas as pd

def bootstrap(invec):
    idx = npr.randint(0, len(invec), len(invec))
    return [invec[i] for i in idx]

def estimatePi0(p, numBoot=100, numLambda=100, maxLambda=0.95):
    p.sort()
    n = len(p)
    lambdas = np.linspace(maxLambda/numLambda,maxLambda,numLambda)
    Wls = np.array([n-np.argmax(p>=l) for l in lambdas])
    pi0s = np.array([Wls[i] / (n * (1 - lambdas[i])) for i in range(numLambda)])
    minPi0 = np.min(pi0s)
    mse = np.zeros(numLambda)
    for boot in range(numBoot):
        pBoot = bootstrap(p)
        pBoot.sort()
        WlsBoot =np.array([n-np.argmax(pBoot>=l) for l in lambdas])
        pi0sBoot =np.array([WlsBoot[i] / (n *(1 - lambdas[i])) for i in range(numLambda)])
        mse = mse + np.square(pi0sBoot-minPi0)
    minIx = np.argmin(mse)
    return pi0s[minIx]

def qvalues(pvalues, pcolname = 'p'):

    m = float(len(pvalues.transpose().values[0].tolist()))
    assert(m>0)
    pvalues = pvalues.sort_values(by = pcolname)
    pi0 = estimatePi0(pvalues.transpose().values[0].tolist())
    num_p, p_sum = 0, 0.0

    qs = pd.DataFrame(columns = [ 'q' ])

    for index, row in pvalues.iterrows():
        p = row[pcolname]
        num_p += 1
        p_sum += p
        q = pi0*p*m/float(num_p)
        qs.loc[index,'q'] = q 


    qs = qs.iloc[::-1]
    old_q=1.0
    for ix in range(len(qs)):
        q = min(old_q,qs.iloc[ix, 0])
        old_q = q
        qs.iloc[ix, 0] = q
    return qs