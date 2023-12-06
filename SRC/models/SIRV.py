import numpy as np
import igraph as ig
from Clustering import *


import sys
sys.path.append('..')
from utilities.IC import *



def init_SIRV(P_dyn, G):
    """
    SIRV is initialized with the following parameters:
    "Dynamic_X" :{
            "func" : "SIRV",
            "S2I"   : 0.2,
            "I2R"   : 0.1,
            "S2V_l" : 0,
            "S2V_h" : 0.2,
            "VA" : 0.1,                         # Vaccines Availability: maximal fraction of the population that can be vaccinated in one time step
            "HEALTH" : {"layer" : 0, "IC" : {"type" : "random", "N_pat_zero" : 1}},
            "BEHAVIOR" : {"layer" : 0, "static" : "True", "IC" :   { "type" : "random_beta", "alpha" : 0.1,
                                                  "homophily" : { "Flag" : "True", "target" : -1, "steps" : 10000, "precision" : 0.01}}                                                                            }
    }
    """
    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    G[hl].vs["S2I"] = P_dyn["S2I"]                                                              # for each node sets  beta
    G[hl].vs["I2R"]   = P_dyn["I2R"]                                                            # for each node sets  gamma
    G[hl].vs["S2V_l"] = P_dyn["S2V_l"]                                                          # for each node sets the low  vacc probability
    G[hl].vs["S2V_h"] = P_dyn["S2V_h"]                                                          # for each node sets the high vacc probability
    G[hl].vs["VA"]    = P_dyn["VA"]                                                             # for each node sets the availability of vaccines
    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    set_continuous_initial_condition(P_dyn["BEHAVIOR"]["IC"], "behavior_status", G[bl])

    if(P_dyn["BEHAVIOR"]["static"]):
        PB = np.array(G[bl].vs["behavior_status"])
        # P(V) = P_l + PB * (P_h - P_l)
        G[hl].vs["chances_to_vaccinate"] = G[hl].vs["S2V_l"][0] + PB * (G[hl].vs["S2V_h"][0] - G[hl].vs["S2V_l"][0])    # if this is static, I calculate it only once, and then I use it in the update function


    rule  = {
        'rule': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'static' : P_dyn["BEHAVIOR"]["static"]
    }
    return rule



def update_SIRV(G, rule):
    #input("updating SIRV")

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
            g_b.vs["chances_to_vaccinate"] = g_h.vs["S2V_l"][0] + PB * (g_h.vs["S2V_h"][0] - g_h.vs["S2V_l"][0])

        S2I = g_h.vs["S2I"][0] # it's the same for everyone, I put it here to save time, faster than accessing the dictionary every time
        I2R = g_h.vs["I2R"][0] # it's the same for everyone, I put it here to save time, faster than accessing the dictionary every time
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(len(g_h.vs)))

        # first, people who wants to vaccinate will try to do so.
        for i,vertex in enumerate(g_h.vs):
            # only a susceptible person can vaccinate
            if vertex["health_status"]==1:
                if(rr[i]<g_h.vs[i]["chances_to_vaccinate"]):
                    vertex['temp_status']=4
        # these are people who are trying to be vaccinated. if more than VA * N people are trying to be vaccinated, only VA * N will be vaccinated, they are picked at random

        indx_vacc = np.where(np.array(g_h.vs["temp_status"])==4)[0]
        if(len(indx_vacc)>g_h.vs["VA"][0]*len(g_h.vs)):
            # if more than VA * N people are trying to be vaccinated, only VA * N will be vaccinated, they are picked at random
            #select a random sample from indx_vacc of len = len(indx_vacc) - VA * N # they will fail to get vaccinated
                                                                                    # and their status will be set back to 1
            L = len(indx_vacc) - int(g_h.vs["VA"][0]*len(g_h.vs))
            for gg in np.random.choice(indx_vacc, L, replace = False):
                g_h.vs[gg]['temp_status']=1
        # now, a normal update of the SIR model
        for i,vertex in enumerate(g_h.vs):
            if vertex["temp_status"]==4:
                vertex["health_status"]=4

        for i,vertex in enumerate(g_h.vs):
            if vertex["health_status"]==1:

                h_n = np.sum(np.array(g_h.vs[g_h.neighbors(i)]["health_status"])==2)
                if(h_n>0):
                    pinf = 1-np.power(1-S2I,h_n)
                    if(rr[i]<pinf):
                        vertex['temp_status']=2
                    else:
                        vertex['temp_status']=1
                else:
                    vertex['temp_status']=1
            elif(vertex["health_status"]==2):
                if(rr[i]<I2R):
                    vertex['temp_status']=3
                else:
                    vertex['temp_status']=2
            elif vertex["health_status"]==3:
                vertex['temp_status']=3
            elif vertex["health_status"]==4:
                vertex['temp_status']=4
            else:
                print("state = ", vertex['temp_status'], vertex["health_status"])
                input('ERROR! This node had an health state that was unexpected!!! in evolve_SIRB_beta')

        # update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"]=vertex['temp_status']

def init_model(update_dictionary, init_single_rule):
    update_dictionary["SIRV"] = update_SIRV
    init_single_rule["SIRV"] = init_SIRV