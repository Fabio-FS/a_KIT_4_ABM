#!/usr/bin/env python
# coding: utf-8

# In[ ]:
# meaningless change to test git2

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



def calc_polarization(g, name = "behavior_status"):
    B = np.array(g.vs[name])
    return polarization(B)

def polarization(B):
    A2= np.ones([len(B),len(B)])- np.diag(np.ones(len(B)))
    pol = 0
    for i in range(len(B)):
        pol = pol + np.sum(A2[i,:].dot(np.power(B[i]-B,2)))



def calc_homophily(g, name = "behavior_status"):    # function called by the save_homophily function in Save_Functions.py
    A = np.array(g.get_adjacency().data)
    B = np.array(g.vs[name])
    H = homophily(B,A)                              # function defined below
    return H

def homophily(B,A):
    H = homophily_non_rescaled(B,A)
    resc = rescale_H(B,A, N_trials = 100)
    H = H + resc
    m2 = 2/np.sum(A)
    H = H*m2
    return H

def homophily_non_rescaled(B,A, metric = "L2"):
    H = 0
    for i in range(len(B)):      
        H = H + np.sum(A[:,i].dot(0.5-np.power(B[i]-B,2)))
    return H

def rescale_H(B,A, N_trials = 1000):
    res = 0
    for i in range(N_trials):
        B_rand = np.random.permutation(B)
        H = homophily_non_rescaled(B_rand, A)
        res = res + H
    res = res/N_trials
    return -res






def metropolis(g, name =  "behavior_status", target = 0, N_steps = 10, return_H_hist = False, dbg = False):

    # initialize the history of all the behaviors, this IS VERY MEMORY INTENSIVE. for each time-step of the metropolis algorithm, we store the behavior of all the nodes.
    if(dbg):
        hist_B = np.zeros([N_steps,len(g.vs)])
    
    Hs, B, A, L, m2 = initialize(g, name, N_steps)  # history of homophily, vector of behaviors, adjacency matrix, number of nodes, 2/sum(A)
    count = 0
    k12 = np.random.choice(np.arange(L), [N_steps,2])

    Hs[0] = homophily(B,A)/m2
    
    target = target/m2

    while(count < N_steps):

        i, j = k12[count,0], k12[count,1]
        B2 = B.copy()
        dh = -calc_dh(i,j,B2,A)
        B2[i], B2[j] = B2[j], B2[i]

        h_attempt = Hs[count]+dh

        
        if(np.abs(Hs[count]-target)>=np.abs(h_attempt-target)):
            Hs[count+1] = h_attempt
            B[i], B[j] = B[j], B[i]
        else:
            Hs[count+1] = Hs[count]
        if(dbg):
            hist_B[count,:] = B
        count = count+1
    Hs = Hs*m2

    

    if (not return_H_hist):
        if(not dbg):
            results = (Hs[-1])
    if (return_H_hist):

        if(not dbg):
            results = (Hs[-1],Hs)
    if (not return_H_hist):

        if(dbg):
            results = (Hs[-1],hist_B, A)
    if (return_H_hist):
        if(dbg):
            results = (Hs[-1],Hs, hist_B, A)
    
    return results
    


def initialize(g, name, N_steps):
    Hs = np.ones(N_steps+1)      # history of homophily
    B = np.array(g.vs[name])
    A = np.array(g.get_adjacency().data)
    L = len(B)
    m2 = 2/np.sum(A)
    return Hs, B, A, L, m2
    

def calc_dh(k,l,B,A):

    B2 = B.copy()
    B2[k], B2[l] = B2[l], B2[k]

    dh = A[:,k].dot(
        np.power(B2-B2[k],2)-np.power(B-B[k],2)) + A[:,l].dot(
        np.power(B2-B2[l],2)-np.power(B-B[l],2))

    return dh*2