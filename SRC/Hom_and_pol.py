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

# This is a Metropolis algorithm, each proposed move switches the behavior of two nodes.
# CERCA tra gli appunti di DI RENZO le note sugli algoritmi per ISING



def calc_polarization(g, name = "behavior_status"):
    B = np.array(g.vs[name])
    return polarization(B)

def polarization(B):
    A2= np.ones([len(B),len(B)])- np.diag(np.ones(len(B)))
    pol = 0
    for i in range(len(B)):
        pol = pol + np.sum(A2[i,:].dot(np.power(B[i]-B,2)))
    return pol/np.sum(A2)



def calc_homophily(g, name = "behavior_status", flag = 0):    # function called by the save_homophily function in Save_Functions.py
    A = np.array(g.get_adjacency().data)
    B = np.array(g.vs[name])
    H = homophily_non_rescaled(B,A)
    if(flag == 0):                              # when I call the function with flag = 0, I want to calculate the homophily and rescale it only if the flag is set to 1 in the graph.
        if(g["rescale_homophily_flag"] == 1):
            H = H + rescale_H(B,A)                      
    elif(flag == 1):                            # when I call the function with flag = 1, I want to calculate the homophily and rescale it.
        H = H + rescale_H(B,A)
    elif(flag == 2):                            # when I call the function with flag = 2, I want to calculate the homophily without rescaling it.
        H = H    
    return H


def homophily_non_rescaled(B,A):
    H = 0
    for i in range(len(B)):      
        H = H + np.sum(A[:,i].dot(0.5-np.power(B[i]-B,2)))
    H = 2*H/np.sum(A)
    return H

def rescale_H(B,A, N_trials = 100):
    L = len(B)
    A2 = np.ones([L,L])-np.diag(np.ones(L))
    dummy_resc = homophily_non_rescaled(B,A2)              # very good approximation of the rescaling factor. It's much faster than the exact calculation.
        
    #res = 0
    #for i in range(N_trials):
    #    B_rand = np.random.permutation(B)
    #    H = homophily_non_rescaled(B_rand, A)
    #    res = res + H
    #res = res/N_trials
    return -dummy_resc






def metropolis(g, name =  "behavior_status", target = 0, N_steps = 10, tollerance = 1e-4, return_H_hist = False, dbg = False, rescale = False):

    # initialize the history of all the behaviors, this IS VERY MEMORY INTENSIVE. for each time-step of the metropolis algorithm, we store the behavior of all the nodes.
    if(dbg):
        hist_B = np.zeros([N_steps,len(g.vs)])
    
    Hs, B, A, L, m2 = initialize(g, name, N_steps)  # history of homophily, vector of behaviors, adjacency matrix, number of nodes, 2/sum(A)
    count = 0
    k12 = np.random.choice(np.arange(L), [N_steps,2])

    Hs[0] = homophily_non_rescaled(B,A)
    #print("Initial homophily: ", Hs[0], "Target: ", target, "m2: ", m2)
    if (rescale):
        resc = rescale_H(B,A)
        Hs[0] = Hs[0] + resc
        g["rescale_homophily_flag"] = 1         # I add a flag to the graph to remember that I rescaled the homophily. It's quite ugly, but it's the only way I found to keep track of it.
                                                        # it is needed in the calc_homophily function.
        
        ### this part is for testing Sven idea. much faster if it works well.
        #print("real_resc = ", resc, "dummy_resc = ", dummy_resc, "difference = ", resc-dummy_resc, "percentage = ", (resc-dummy_resc)/resc*100, "%")


    else:
        g["rescale_homophily_flag"] = 0

    #g["homophily0"] = Hs[0]
#    print("Initial homophily: ", Hs[0])
#    print("Target: ", target)
#    target = target

    while(count < N_steps):
        #if(count%10000 == 0):
            #print("Step: ", count)
        #error = np.abs(Hs[count]-target)
        #print("Error: ", error)
        if(tollerance < np.abs(Hs[count]-target)):
            #print("try")
            i, j = k12[count,0], k12[count,1]
            B2 = B.copy()
            dh = -calc_dh(i,j,B2,A)*m2
            B2[i], B2[j] = B2[j], B2[i]

            h_attempt = Hs[count]+dh

            
            if(np.abs(Hs[count]-target)>=np.abs(h_attempt-target)):
                #print("Accepted move: ", count, ", H moved of ", h_attempt-Hs[count])  
                Hs[count+1] = h_attempt
                B[i], B[j] = B[j], B[i]
            else:
                #print("Rejected move: ", count, "  Error: ", np.abs(Hs[count]-target), "  Attempt: ", np.abs(h_attempt-target))
                Hs[count+1] = Hs[count]
            if(dbg):
                hist_B[count,:] = B
            count = count+1
        else:
            Hs[count+1] = Hs[count]
            count = count +1
    Hs = Hs
    g.vs[name] = B

    

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
    
    # I print the final homophily, and I check if the function homophily works correctly:
    #print("Target: ", target, "Final homophily: ", Hs[-1], " Delta: ", Hs[-1]-target)


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