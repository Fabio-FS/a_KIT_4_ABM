import numpy as np
import igraph as ig
from Clustering import *
from Metrics_of_graph import *
import os


def set_disease_initial_condition(IC, name, g):
    if(IC["type"] == "random"):
        g.vs[name]=1
        # select N_pat_zero nodes at random and set them to 2:
        for gg in rn.sample(range(len(g.vs)), IC["N_pat_zero"]):

            g.vs[gg][name] = 2
    elif(IC["type"] == "central_geometric"):
        g.vs[name]=1
        # select N_pat_zero nodes at random and set them to 2:
        gg = int(np.floor(len(g.vs)/2))
        g.vs[gg][name] = 2
    elif(IC["type"] == "vector"):
        g.vs[name] = IC["values"]
    else:
        print("INITALIZATION RULE: " + IC["type"] + " NOT IMPLEMENTED YET! NUUUUUU")



# implementation of the initial conditions. Used for continuous dynamics

def set_continuous_initial_condition(IC, name, g):
    

    if(IC["type"] == "random_uniform"):
        g.vs[name]=np.random.uniform(IC["Low"],IC["High"],len(g.vs))
    elif(IC["type"] == "random_beta"):
        # if alpha is not defined, sets it to be equal to beta, and viceversa
        if("alpha" not in IC):
            a = IC["beta"]
            b = IC["beta"]
        elif("beta" not in IC):
            a = IC["alpha"]
            b = IC["alpha"]
        else:
            a = IC["alpha"]
            b = IC["beta"]
        #print("getting here!")
        VECTOR_INIT = np.random.beta(a,b,len(g.vs))
        g.vs[name]= VECTOR_INIT
    elif(IC["type"] == "vector"):
        g.vs[name] = IC["values"]
    else:
        print("INITALIZATION RULE: " + IC["type"] + " NOT IMPLEMENTED YET! NUUUUUU")
    if(IC["homophily"]["Flag"] == "True"):
        # case 1 debug = False, return_H_hist = False
        if(IC["homophily"]["debug"] == False and IC["homophily"]["return_H_hist"] == False):
            metropolis_slow(g, name =  name, target = IC["homophily"]["target"], N_steps = IC["homophily"]["steps"], precision = 0.01, return_H_hist = False,debug = False)
            #don't save anything
        # case 2 debug = False, return_H_hist = True
        elif(IC["homophily"]["debug"] == False and IC["homophily"]["return_H_hist"] == True):
            H_hist = metropolis_slow(g, name =  name, target = IC["homophily"]["target"], N_steps = IC["homophily"]["steps"], precision = 0.01, return_H_hist = True,debug = False)
            # save the history of homophily
            
            if os.path.exists("H_hist.csv"):
                os.remove("H_hist.csv")
            np.savetxt("H_hist.csv", H_hist, delimiter=",")
        elif(IC["homophily"]["debug"] == True and IC["homophily"]["return_H_hist"] == False):
            pass # this case is not implemented
        else:
            H_hist, all_B, A = metropolis_slow(g, name =  name, target = IC["homophily"]["target"], N_steps = IC["homophily"]["steps"], precision = 0.01, return_H_hist = True,debug = True)
            # if "H_hist.csv" exists, delete it
            # imports the libraries needed to os.remove
            if os.path.exists("H_hist.csv"):
                os.remove("H_hist.csv")
            # save the history of homophily
            np.savetxt("H_hist.csv", H_hist, delimiter=",")
            # save the history of the behavior as TxM matrix:

            if os.path.exists("all_B.csv"):
                os.remove("all_B.csv")
            

            # decrease the number of all_B to be written on file.
            indx = calc_indxs(all_B)
        
            all_B2 = all_B[indx,:]

            np.savetxt("all_B.csv", all_B2, delimiter=",")

            # save the adjacency matrix
            if os.path.exists("ADJ.csv"):
                os.remove("ADJ.csv")
            ig.Graph.write_adjacency(g, "ADJ.csv")


def calc_indxs(all_B):
    indx = []
    for i in range(int(np.log10(len(all_B)))):
        # append np.power(10,i) to indx
        indx.append(np.power(10,i))
        indx.append(2*np.power(10,i))
        indx.append(4*np.power(10,i))
    indx.append(len(all_B)-1)
    return indx