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

#                                                                                                             
#  █████     ████    ██         ██     █████     ████    ██████     ██     ██████    ████     ████    ██  ██  
#  ██  ██   ██  ██   ██        ████    ██  ██     ██         ██    ████      ██       ██     ██  ██   ███ ██  
#  ██  ██  ██    ██  ██       ██  ██   ██  ██     ██        ██    ██  ██     ██       ██    ██    ██  ██████  
#  █████   ██    ██  ██       ██████   █████      ██       ██     ██████     ██       ██    ██    ██  ██████  
#  ██      ██    ██  ██       ██  ██   ████       ██      ██      ██  ██     ██       ██    ██    ██  ██ ███  
#  ██       ██  ██   ██       ██  ██   ██ ██      ██     ██       ██  ██     ██       ██     ██  ██   ██  ██  
#  ██        ████    ██████   ██  ██   ██  ██    ████    ██████   ██  ██     ██      ████     ████    ██  ██  
#

def calc_polarization(g, attribute = "behavior_status"):
    B = np.array(g.vs[attribute])
    return polarization(B)

def polarization(B):
    A2= np.ones([len(B),len(B)])- np.diag(np.ones(len(B)))
    pol = 0
    for i in range(len(B)):
        pol = pol + np.sum(A2[i,:].dot(np.power(B[i]-B,2)))
    return pol/np.sum(A2)

#                                                                                  
#  ██  ██    ████    ██   ██   ████    █████    ██  ██    ████    ██       ██  ██  
#  ██  ██   ██  ██   ███ ███  ██  ██   ██  ██   ██  ██     ██     ██       ██  ██  
#  ██  ██  ██    ██  ███████ ██    ██  ██  ██   ██  ██     ██     ██       ██  ██  
#  ██████  ██    ██  ██ █ ██ ██    ██  █████    ██████     ██     ██        ████   
#  ██  ██  ██    ██  ██   ██ ██    ██  ██       ██  ██     ██     ██         ██    
#  ██  ██   ██  ██   ██   ██  ██  ██   ██       ██  ██     ██     ██         ██    
#  ██  ██    ████    ██   ██   ████    ██       ██  ██    ████    ██████     ██    
#     

def calc_homophily(g, attribute = "behavior_status", recenter_flag = False, quantity_scaling = 1, is_category = 0):    # function called by the save_homophily function in Save_Functions.py
    #category is an integer misused as a boolean.
    #   It describes if homophily should be evaluated in terms of metric distance or categories

    #quantity_scaling is a scaling factor for the quantity

    #flag describes if the homophily value should be recentered.
    #   When all nodes have the same value, no amount of rewiring can change the homophily.
    #   Depending on the usecase, this might be desirable to interpret as "no homophily (above expected)"
    #   If every node but one has the same value. With the odd one being isolate, this would be a higher homophily than was expected just based on the feature vector.

    A = np.array(g.get_adjacency().data)
    if type(attribute) == np.ndarray:
        B = attribute
    else:
        B = np.array(g.vs[attribute])

    hom_n_rc = [homophily_non_recentered,cat_homophily_non_recentered][is_category]
    rc_H = [recenter_H,cat_recenter_H][is_category]

    H = hom_n_rc(B,A, quantity_scaling = quantity_scaling)                    
    if recenter_flag == True :                            # when I call the function with flag = 1, I want to calculate the homophily and recenter it.
        H = H + rc_H(B,A, quantity_scaling = quantity_scaling)
    else:                            # when I call the function with flag = 2, I want to calculate the homophily without recentering it.
        H = H    
    return H

def homophily_non_recentered(B,A,quantity_scaling=1):
    #A is the adjacency matrix, B the vector of quantity of interest
    #If the quantity is \element [0,0.2], quantity_scaling should be 5
    #if it is \element [3,7], quantity_scaling should be 0.25

    C = np.tile(B,(B.shape[0],1)).T
    #C is a tiling of the attribute vector, but  transposed
    # ex: B = [0,1,2] -> C = [[0 0 0] [1 1 1] [2 2 2]]
    
    H = np.sum(  A * (0.5 - (quantity_scaling*(C-B))**2)   )

    return 2*H/np.sum(A)

def cat_homophily_non_recentered(B,A,quantity_scaling = 1):
    #calculating a categorical homophily, where all values have the same distance to each other
    #A is the adjacency matrix, B the vector of quantity of interest
    #I don't need quantity_scaling here, but I have it here to be interchangable with the other homophily function

    C = np.tile(B,(B.shape[0],1)).T
    #C is a tiling of the attribute vector, but  transposed
    # ex: B = [0,1,2] -> C = [[0 0 0] [1 1 1] [2 2 2]]

    H = np.sum( A * ((C==B) - 1/2) )
    #(C==B)  is a matrix with a 1 at index i,j  if  Bi == Bj;  otherwise all zeros
    #from that matrix, you subtract half.  so now C_ij = 0.5 if Bi == Bj;   otherwise -0.5
    #multiply with the adjacency matrix to filter out only connected nodes:    A* ((C==B) - 1/2)
    #and sum all values
    #each connection of same agents adds 1/2, each connection of different agents subtracts 1/2
    
    return 2*H/np.sum(A)

def recenter_H(B,A, quantity_scaling = 1):
    L = len(B)
    A2 = np.ones([L,L])-np.diag(np.ones(L))    #adjacency matrix of complete graph
    dummy_resc = homophily_non_recentered(B,A2, quantity_scaling = quantity_scaling)              # very good approximation of the rescaling factor. It's much faster than the exact calculation.

    return -dummy_resc

def cat_recenter_H(B,A, quantity_scaling = 1):
    L = len(B)
    A2 = np.ones([L,L])-np.diag(np.ones(L))    #adjacency matrix of complete graph
    recenter_value = cat_homophily_non_recentered(B,A2)              # very good approximation of the recenter value. It's much faster than the exact calculation.

    return -recenter_value



def calculate_homophily_difference(k,l,B,A):

    #calculate difference in homophily if two nodes   k and l   were switched.
    #B: the attribute vector before the switch
    #A: the adjacency matrix

    if A[k,l] == 0:
        #not neighbors - very likely in sparsely connected graphs
        #then there is no need to make copies
        dh = ( A[:,k].dot( (B-B[l])**2 - (B-B[k])**2 ) +
               A[:,l].dot( (B-B[k])**2 - (B-B[l])**2 )  )
        return 2*dh
    
    #if the two nodes are neighbors, the change in homophily needs to be calculated differently.
    #then the above formula overestimates homophily because one of k's neighbors IS l.
    #but I'm calculating how similar l would be to k's neighbors, one of which is l themselves
    #that's why in this case, I just copy the whole vector and make the switch.

    B2 = B.copy()
    B2[k], B2[l] = B2[l], B2[k]

    dh = ( A[:,k].dot( (B2-B2[k])**2 - (B-B[k])**2 ) +
           A[:,l].dot( (B2-B2[l])**2 - (B-B[l])**2 )  )

    return 2*dh

def cat_calculate_homophily_difference(k,l,B,A):

    #calculate difference in homophily if two nodes   k and l   were switched.
    #B: the attribute vector before the switch;
    #   the values of B are categories, not continuous
    #A: the adjacency matrix

    if B[k] == B[l]:
        return 0
    
    if A[k,l] == 0:
        #not neighbors - very likely in sparsely connected graphs
        #then there is no need to make copies
        dh = ( A[:,k].dot( (B==B[l]).astype(int) - (B == B[k])) + 
               A[:,l].dot( (B==B[k]).astype(int) - (B == B[l]))    )
        return 2*dh

    #if the two nodes are neighbors, the change in homophily needs to be calculated differently.
    #then the above formula overestimates homophily because one of k's neighbors IS l.
    #but I'm calculating how similar l would be to k's neighbors, one of which is l themselves
    #that's why in this case, I just copy the whole vector and make the switch.
    
    B2 = B.copy()
    B2[k], B2[l] = B2[l], B2[k]

    dh = ( A[:,k].dot( (B2==B[l]).astype(int) - (B == B[k])) + 
           A[:,l].dot( (B2==B[k]).astype(int) - (B == B[l]))    )

    return 2*dh

#                                                                                           
#  ██   ██  ██████   ██████   █████     ████    █████     ████    ██        ████     ████   
#  ███ ███  ██         ██     ██  ██   ██  ██   ██  ██   ██  ██   ██         ██     ██  ██  
#  ███████  ██         ██     ██  ██  ██    ██  ██  ██  ██    ██  ██         ██     ██      
#  ██ █ ██  ████       ██     █████   ██    ██  █████   ██    ██  ██         ██      ████   
#  ██   ██  ██         ██     ████    ██    ██  ██      ██    ██  ██         ██         ██  
#  ██   ██  ██         ██     ██ ██    ██  ██   ██       ██  ██   ██         ██     ██  ██  
#  ██   ██  ██████     ██     ██  ██    ████    ██        ████    ██████    ████     ████   
#

def metropolis(g, attribute_vector,
               hom_target = 0, tolerance = 1e-4,
               recenter = False,
               N_steps = 10,
               is_category = 0, temperature = 0.01,  quantity_scaling = 1,              
               return_H_hist = False, dbg = False,               
               global_var = None
               ):
    # This is a Metropolis algorithm, each proposed move switches the behavior of two nodes.
    #
    # from the graph    g   , we calculate the adjacency matrix A   [shape (N,N)]
    #    attribute     is the vector, that should be redistributed amongst the nodes [shape N]
    #
    #all other values are optional.
    #hom_target:   the target value for homophily
    #tolerance:   how close we need to get to the target - break condition
    #recenter:   factor in expected homophiily given the attribute_vector
    #N_steps:   the maximum number of switching steps
    #is_category:    0 or 1. when True, two different values always have a distance of 1
    #temperature:    how noisy is the acceptance of proposed switches
    #quantity_scaling:   scales the vector. Differences between two values should be never higher than 1.
    #return_H_hist:    whether to write all homophily values to a csv
    #dbg:    whether to return the whole history of attribute_vectors
    #global_var:    I drag this through the whole process to keep track of stuff


    # initialize the history of all the behaviors, this IS VERY MEMORY INTENSIVE. for each time-step of the metropolis algorithm, we store the behavior of all the nodes.
    if(dbg):
        hist_B = np.zeros([N_steps+1,len(g.vs)])
    
    Hs = np.ones(N_steps+1)      # history of homophily
    
    B = attribute_vector
    A = np.array(g.get_adjacency().data)  #adjacency matrix
    edge_weight = 2/np.sum(A)

    count = 0
    random_indcs = np.random.choice(np.arange(g.vcount()), [N_steps,2])

    hom_n_rc = [homophily_non_recentered,cat_homophily_non_recentered][is_category]
    rc_H = [recenter_H,cat_recenter_H][is_category]
    calc_dh = [calculate_homophily_difference, cat_calculate_homophily_difference][is_category]
    Hs[0] = hom_n_rc(B,A, quantity_scaling = quantity_scaling)

    if recenter:
        Hs[0] += rc_H(B,A, quantity_scaling = quantity_scaling)

    break_count = N_steps

    for count in range(N_steps):
        if(tolerance < np.abs(Hs[count]-hom_target)):
            #if tolerance isn't yet reached

            i, j = random_indcs[count,0], random_indcs[count,1]
            dh = calc_dh(i,j,B,A)*edge_weight

            h_attempt = Hs[count]+dh

            hom_difference = np.abs(Hs[count]-hom_target) - np.abs(h_attempt-hom_target)
            #how close the algorithm got in the previous step minus how close we are this time
            #if hom_difference > 0:    new value is closer
            #if hom_difference < 0:    old value is closer

            if temperature == 0:
                if hom_difference >= 0 :
                    Hs[count+1] = h_attempt
                    B[i], B[j] = B[j], B[i]
                else:
                    Hs[count+1] = Hs[count]
            else:
                if np.random.random() < 1 / ( 1   +   np.exp( - (1/temperature) * hom_difference)    ) :
                    Hs[count+1] = h_attempt
                    B[i], B[j] = B[j], B[i]
                else:
                    Hs[count+1] = Hs[count]

            if(dbg):
                hist_B[count+1,:] = B
        else:
            Hs[count:] = Hs[count]
            break_count = count
            break

    if not global_var is None:
        setattr(global_var,"metropolis_steps",break_count)
    

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

    setattr(global_var,"homophily_recentered",float(Hs[-1]))   if recenter    else     setattr(global_var,"homophily",float(Hs[-1]))

    return B, results
    
