import numpy as np
import igraph as ig
from Clustering import *


import sys
#sys.path.append('..')
from utilities.IC import *


def init_SIRB_beta(P_dyn, G):
    """
    SIRB_beta is initialized with the following parameters:
    "Dynamic_X" :{
            "func" : "SIRB_beta",
            "S2I_l" : 0,
            "S2I_h" : 0.2,
            "I2R"   : 0.1,
            "HEALTH" : {"layer" : 0, "IC" : {"type" : "random", "N_pat_zero" : 1}},
            "BEHAVIOR" : {"layer" : 0, "IC" :   { "type" : "random_beta", "alpha" : 0.1,
                                                  "homophily" : { "Flag" : "True", "target" : -1, "steps" : 10000, "precision" : 0.01}}                                                                            }
    }
    """
    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    G[hl].vs["S2I_h"] = P_dyn["S2I_h"]                                                          # for each node sets the high beta
    G[hl].vs["S2I_l"] = P_dyn["S2I_l"]                                                          # for each node sets the low beta
    G[hl].vs["I2R"]   = P_dyn["I2R"]                                                            # for each node sets the gamma

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    set_continuous_initial_condition(P_dyn["BEHAVIOR"]["IC"], "behavior_status", G[bl])

    if(P_dyn["BEHAVIOR"]["static"]):
        PB = np.array(G[bl].vs["behavior_status"])
        Low  = G[hl].vs["S2I_l"][0]
        High = G[hl].vs["S2I_h"][0]
        G[hl].vs["updated_beta"] = Low + PB*(High-Low)    # if this is static, I calculate it only once, and then I use it in the update function

    #for i in range(len(G[hl].vs)):
        #print("node ", i, " beta_i = ", np.round(G[hl].vs[i]["updated_beta"],2), " and behavior", np.round(G[bl].vs[i]["behavior_status"],2) )

    rule  = {
        'rule': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'static' : P_dyn["BEHAVIOR"]["static"]
    }
    return rule



def update_SIRB_beta(G, rule):
    """this simulates the classical SIR model, where the probability of being infected depends ONLY on the behavior of the susceptible individuals
    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer]
    """
    g_h = G[rule["hl"]]
    g_b = G[rule["bl"]]
    #name = g_h.vs["health_status"][0]
    # check if there are infected nodes, if not, skip the update
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)

    if(N_infected>0):
        if(rule["static"]):
            pass
        else:
            # create an array with the susceptibility of each node
            PB = np.array(g_b.vs["behavior_status"])
            G[g_h].vs["updated_beta"] = g_h.vs["S2I_l"][0] + PB*(g_h.vs["S2I_h"][0]-g_h.vs["S2I_l"][0])

        I2R = g_h.vs["I2R"][0] # it's the same for everyone, I put it here to save time, faster than accessing the dictionary every time
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(len(g_h.vs)))
        for i,vertex in enumerate(g_h.vs):
            # the updated health status is not overwritten immediately into health_status, but it is stored in temp_status to avoid updating the health status of a
            # subsequent node with an updated value. this keeps the update syncronous.
            
            # for each node i, if it is susceptible, check if it gets infected:
            if vertex["health_status"]==1:
                h_n = np.sum(np.array(g_h.vs[g_h.neighbors(i)]["health_status"])==2)
                #print("I have ", h_n, " infected neighbors")
                if(rr[i]<1-np.power(1-g_h.vs[i]["updated_beta"],h_n)):
                    vertex['temp_status']=2
                else:
                    vertex['temp_status']=1

            elif(vertex["health_status"]==2):
                if(rr[i]<I2R):
                    vertex['temp_status']=3
                else:
                    vertex['temp_status']=2
            elif vertex["health_status"]==3:
                vertex['temp_status']=3
            else:
                input('ERROR! This node had an health state that was unexpected!!! in evolve_SIRB_beta')

        # update the health status of the nodes
        for vertex in g_h.vs:
            #print("old health status: ", vertex["health_status"], " new health status: ", vertex['temp_status'])
            vertex["health_status"]=vertex['temp_status']


def init_model(update_dictionary, init_single_rule):
    update_dictionary["SIRB_beta"] = update_SIRB_beta
    init_single_rule["SIRB_beta"] = init_SIRB_beta

