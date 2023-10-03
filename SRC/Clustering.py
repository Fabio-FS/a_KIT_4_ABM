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


def Metropolis_single_switch(g, name, target = 0, N_steps = 1000, precision = 0.01):

    
    A = np.array(g.get_adjacency().data)
    B = np.array(g.vs[name])
    
    BO = B
    NL = 2/np.sum(A)
    L = len(g.vs)
    
    Hs=np.zeros(N_steps+1)

    OH = exact_H(A,B)
    
    
    count = 0
    k12 = np.random.choice(np.arange(L), [N_steps,2])
    
    A = A * NL 
    
    while(count < N_steps):
        Hs[count] = OH
        
        i = k12[count,0]
        j = k12[count,1]
        
        dp = update_dp(A,B,i,j)
        NH =  OH + dp
        DH = np.abs(OH-target)-np.abs(NH-target)
    
        if(DH>=0):
            OH = NH
        
        else:
            B[i], B[j] = B[j], B[i]
        count = count+1
        if (np.abs(OH-target)<precision):
            #print("quit for shortcut")
            for t in range(count,N_steps):
                Hs[t] = OH
            g.vs["behavior_status"] = B
            return Hs, count
    g.vs["behavior_status"] = B
    
    Hs[count] = OH        
    
    return Hs




def update_dp(A,B,i,j):
    dp = 0
    
    Vi = np.abs(B-B[i])
    Vj = np.abs(B-B[j])
    dp = -A[:,i].dot(Vi) - A[:,j].dot(Vj)
    
    B[i], B[j] = B[j], B[i]
    
    Vi = np.abs(B-B[i])
    Vj = np.abs(B-B[j])
    dp = dp + A[:,i].dot(Vi) + A[:,j].dot(Vj)
        
    dp = -dp
        
    return dp

def exact_pol(A,B):
    
    NL = 2/np.sum(A)
    Op = 0
    for i in range(len(B)):
        Op = Op - np.sum(A[:,i].dot(np.abs(B[i]-B)))
    
    Op = 1 + ( Op/2)*NL
    H  = Op
    
    return H

def exact_H(A,B):
    A2= np.ones([len(B),len(B)])- np.diag(np.ones(len(B)))

    p = exact_pol(A,B)
    Ep = exact_pol(A2,B)
    
    H = p - Ep
    
    
    return H

def calc_H(g, name = "behavior_status"):
    A = np.array(g.get_adjacency().data)
    B = np.array(g.vs[name])
    H = exact_H(A,B)
    return H


def cluster_behavior(g_beh,CLUSTER, debug = False):
    
    if((CLUSTER.algorithm == "Metropolis_switch") & (debug == False)):
        Metropolis_single_switch(g_beh,CLUSTER)
    elif((CLUSTER.algorithm == "Metropolis_switch") & (debug == True)):
        H_hist = Metropolis_single_switch(g_beh, CLUSTER, debug = True)
        return H_hist
    else:
        print("WARNING! CLUSTERING ALGORITHM " + algorithm + " NOT IMPLEMENTED YET!")
        print("proceeding with Metropolis_single_switch algorithm (only one I have)")
        Metropolis_single_switch(g_beh, CLUSTER)
        
    
    
    
    



