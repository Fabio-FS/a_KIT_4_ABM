#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import igraph as ig
import random as rn
import pathlib
#from numba import jit
import warnings
from scipy.stats import beta

# This is a Metropolis type of algorithm, where the proposed move is to switch the behavior of two nodes.
# Consider, in the future, of implementing also a metropolis move where the behavior of only ONE node is considered.
# CERCA tra gli appunti di DI RENZO le note sugli algoritmi per ISING


def homophily_sq(B,A, metric = "L2"):

    H = 0
    for i in range(len(B)):
        if (metric == "L2"):        
            H = H + np.sum(A[:,i].dot(0.5-np.power(B[i]-B,2)))
        if (metric == "L1"):
            H = H + np.sum(A[:,i].dot(0.5-np.abs(B[i]-B)))

    return H


def rescale_H(B,A, metric = "L2", N_trials = 1000):
    res = 0
    for i in range(N_trials):
        B_rand = np.random.permutation(B)
        H = homophily_sq(B_rand, A, metric = metric)
        res = res + H
    res = res/N_trials
    return -res


def initialize(g, name, N_steps):
    Hs = np.ones(N_steps+1)      # history of homophily
    B = np.array(g.vs[name])
    A = np.array(g.get_adjacency().data)
    L = len(B)
    m2 = 2/np.sum(A)

    return Hs, B, A, L, m2


def metropolis_slow(g, name =  "behavior_status", target = 0, N_steps = 10, precision = 0.01, return_H_hist = False, debug = False):
    
    print("I will do ", N_steps, " steps")


    if(debug == True):
        hist_B = np.zeros([N_steps,len(g.vs)])
    
    Hs, B, A, L, m2 = initialize(g, name, N_steps)  # history of homophily, vector of behaviors, adjacency matrix, number of nodes, 2/sum(A)

    count = 0
    k12 = np.random.choice(np.arange(L), [N_steps,2])


    Hs[0] = homophily_sq(B,A)

    resc =  rescale_H(B,A)*m2
    print("rescaling factor: ", resc)
    target = target/m2 - resc

    while(count < N_steps):

        if(count%1000 == 0):
            print(count)

        i, j = k12[count,0], k12[count,1]
        B2 = B.copy()
        dh = -calc_dh(i,j,B2,A)
        B2[i], B2[j] = B2[j], B2[i]

        #H1 = homophily_sq(B,A)
        #H2 = homophily_sq(B2,A)
        #dh_0 = H2-H1
        #print("Delta dh: ", dh-dh_0)

        h_attempt = Hs[count]+dh

        
        if(np.abs(Hs[count]-target)>=np.abs(h_attempt-target)):
        #    print("accepted move, H: ", Hs[count], " H_attempt: ", h_attempt, " target: ", target)
            Hs[count+1] = h_attempt
            B[i], B[j] = B[j], B[i]
        else:
            Hs[count+1] = Hs[count]
        if(debug == True):
            hist_B[count,:] = B
        count = count+1
        
    Hs = Hs*m2 - resc
    print("final H: ", Hs[count])
    #print("exact current H:", homophily_sq(B,A)*m2 + resc)
    
    if return_H_hist:
        if debug:
            return Hs, hist_B, A
        else:
            return Hs
    else:
        if debug:
            return hist_B, A
    
    

def calc_dh(k,l,B,A):

    B2 = B.copy()
    B2[k], B2[l] = B2[l], B2[k]

    dh = A[:,k].dot(
        np.power(B2-B2[k],2)-np.power(B-B[k],2)) + A[:,l].dot(
        np.power(B2-B2[l],2)-np.power(B-B[l],2))

    return dh*2