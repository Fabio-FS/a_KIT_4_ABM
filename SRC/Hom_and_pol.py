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



def calc_homophily(g, attribute = "behavior_status", flag = 0, qs = 1, is_category = 0):    # function called by the save_homophily function in Save_Functions.py
    #category is an integer misused as a boolean.
    #   It describes if homophily should be evaluated in terms of metric distance or categories
    #qs is a scaling factor for the quantity
    #flag describes if the homophily value should be recentered.
    #   When all nodes have the same value, no amount of rewiring can change the homophily.
    #   Depending on the usecase, this might be desirable to interpret as "no homophily (above expected)"
    #   If every node but one has the same value. With the odd one being isolate, this would be a higher homophily than was expected just based on the feature vector.

    A = np.array(g.get_adjacency().data)
    B = np.array(g.vs[attribute])

    hom_n_rc = [homophily_non_recentered,cat_homophily_non_recentered][is_category]
    rc_H = [recenter_H,cat_recenter_H][is_category]

    H = hom_n_rc(B,A, qs = qs)
    if(flag == 0):                              # when I call the function with flag = 0, I want to calculate the homophily and recenter it only if the flag is set to 1 in the graph.
        if(g["recenter_homophily_flag"] == 1):
            H = H + rc_H(B,A, qs = qs)                      
    elif(flag == 1):                            # when I call the function with flag = 1, I want to calculate the homophily and recenter it.
        H = H + rc_H(B,A, qs = qs)
    elif(flag == 2):                            # when I call the function with flag = 2, I want to calculate the homophily without recentering it.
        H = H    
    return H

def homophily_non_recentered(B,A,qs=1):
    #A is the adjacency matrix, B the vector of quantity of interest
    #If the quantity is \element [0,0.2], qs should be 5
    #if it is \element [3,7], qs should be 0.25

    C = np.tile(B,(B.shape[0],1)).T
    #C is a matrix of the attribute, each row corresponds to one node
    H = np.sum(np.sum(A.T*(0.5 - (qs*(C-B))**2),1))
    return 2*H/np.sum(A)

def cat_homophily_non_recentered(B,A,qs = 1):
    #calculating a categorical homophily, where all values have the same distance to each other
    #A is the adjacency matrix, B the vector of quantity of interest
    #I don't need qs here, but I have it here to be interchangable with the other homophily function

    C = np.tile(B,(B.shape[0],1)).T
    #C is a matrix of the attribute, each row corresponds to one node
    H = np.sum(np.sum(A.T*(0.5 - (1-(C==B))**2),1))
    
    return 2*H/np.sum(A)

def recenter_H(B,A, qs = 1):
    L = len(B)
    A2 = np.ones([L,L])-np.diag(np.ones(L))    #adjacency matrix of complete graph
    dummy_resc = homophily_non_recentered(B,A2, qs = qs)              # very good approximation of the rescaling factor. It's much faster than the exact calculation.

    return -dummy_resc

def cat_recenter_H(B,A, qs = 1):
    L = len(B)
    A2 = np.ones([L,L])-np.diag(np.ones(L))    #adjacency matrix of complete graph
    dummy_resc = cat_homophily_non_recentered(B,A2)              # very good approximation of the rescaling factor. It's much faster than the exact calculation.

    return -dummy_resc




def metropolis(g, attribute =  "behavior_status",
               hom_target = 0, recenter = False,
               N_steps = 10, tolerance = 1e-4,
               return_H_hist = False, dbg = False,
               is_category = 0, temperature = 0.01
               ):
    #N_steps: the maximum number of switching steps

    # initialize the history of all the behaviors, this IS VERY MEMORY INTENSIVE. for each time-step of the metropolis algorithm, we store the behavior of all the nodes.
    if(dbg):
        hist_B = np.zeros([N_steps,len(g.vs)])
    
    Hs, B, A, L, m2 = initialize(g, attribute, N_steps)  # history of homophily, vector of behaviors, adjacency matrix, number of nodes, 2/sum(A)
    count = 0
    k12 = np.random.choice(np.arange(L), [N_steps,2])

    hom_n_rc = [homophily_non_recentered,cat_homophily_non_recentered][is_category]
    rc_H = [recenter_H,cat_recenter_H][is_category]
    Hs[0] = hom_n_rc(B,A)
    #HS[0] = calc_homophily()

    #print("Initial homophily: ", Hs[0], "Target: ", target, "m2: ", m2)
    if (recenter):
        Hs[0] = Hs[0] + rc_H(B,A)
        g["recenter_homophily_flag"] = 1         # I add a flag to the graph to remember that I recentered the homophily. It's quite ugly, but it's the only way I found to keep track of it.
                                                        # it is needed in the calc_homophily function.
        
        ### this part is for testing Sven idea. much faster if it works well.
        #print("real_resc = ", resc, "dummy_resc = ", dummy_resc, "difference = ", resc-dummy_resc, "percentage = ", (resc-dummy_resc)/resc*100, "%")
# test

    else:
        g["recenter_homophily_flag"] = 0

    #g["homophily0"] = Hs[0]
#    print("Initial homophily: ", Hs[0])
#    print("Target: ", target)
#    target = target

    while(count < N_steps):
        #if(count%10000 == 0):
            #print("Step: ", count)
        #error = np.abs(Hs[count]-target)
        #print("Error: ", error)
        if(tolerance < np.abs(Hs[count]-hom_target)):
            #print("try")
            i, j = k12[count,0], k12[count,1]
            B2 = B.copy()
            dh = -calc_dh(i,j,B2,A)*m2
            B2[i], B2[j] = B2[j], B2[i]

            h_attempt = Hs[count]+dh

            hom_difference = np.abs(Hs[count]-hom_target) - np.abs(h_attempt-hom_target)
            #how close the algorithm got in the previous step minus how close we are this time

            if temperature == 0:
                if hom_difference >= 0 :
                    #print("Accepted move: ", count, ", H moved of ", h_attempt-Hs[count])
                    Hs[count+1] = h_attempt
                    B[i], B[j] = B[j], B[i]
                else:
                    #print("Rejected move: ", count, "  Error: ", np.abs(Hs[count]-target), "  Attempt: ", np.abs(h_attempt-target))
                    Hs[count+1] = Hs[count]
            else:
                if np.random.random() < 1 / ( 1   +   np.exp( - (1/temperature) * hom_difference)    ) :
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
    g.vs[attribute] = B

    

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