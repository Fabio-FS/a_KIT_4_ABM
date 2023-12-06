import numpy as np
import igraph as ig
from Clustering import *


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
        # uses my Metropolis algorithm to initialize the homophily
        Hom_hist = Metropolis_single_switch(g, name, target = IC["homophily"]["target"], N_steps = IC["homophily"]["steps"], precision = IC["homophily"]["precision"])