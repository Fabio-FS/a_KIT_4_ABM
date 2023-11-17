import numpy as np
import igraph as ig
from Clustering import *

import time


# implementation of the initialization functions





#  ████████  ██████   ██████  ██      ██████   ██████  ██   ██ 
#     ██    ██    ██ ██    ██ ██      ██   ██ ██    ██  ██ ██  
#     ██    ██    ██ ██    ██ ██      ██████  ██    ██   ███   
#     ██    ██    ██ ██    ██ ██      ██   ██ ██    ██  ██ ██  
#     ██     ██████   ██████  ███████ ██████   ██████  ██   ██ 
#                                                              



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



#   ▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄ 
#  ▐░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
#  ▐░█▀▀▀▀▀▀▀█░▌▀▀▀▀█░█▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀▀▀ 
#  ▐░▌       ▐░▌    ▐░▌     ▐░▌          ▐░▌          ▐░▌       ▐░▌▐░▌          ▐░▌          
#  ▐░▌       ▐░▌    ▐░▌     ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░█▄▄▄▄▄▄▄▄▄ 
#  ▐░▌       ▐░▌    ▐░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
#  ▐░▌       ▐░▌    ▐░▌      ▀▀▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀▀▀ ▐░█▀▀▀▀▀▀▀█░▌ ▀▀▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀▀▀ 
#  ▐░▌       ▐░▌    ▐░▌               ▐░▌▐░▌          ▐░▌       ▐░▌          ▐░▌▐░▌          
#  ▐░█▄▄▄▄▄▄▄█░▌▄▄▄▄█░█▄▄▄▄  ▄▄▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ ▐░▌       ▐░▌ ▄▄▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄▄▄ 
#  ▐░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌▐░▌       ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌
#   ▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀         ▀  ▀▀▀▀▀▀▀▀▀▀▀  ▀▀▀▀▀▀▀▀▀▀▀ 
#                                                                                            







#    #####  ### ######                                   
#   #     #  #  #     #       #####  ###### #####   ##   
#   #        #  #     #       #    # #        #    #  #  
#    #####   #  ######  ##### #####  #####    #   #    # 
#         #  #  #   #         #    # #        #   ###### 
#   #     #  #  #    #        #    # #        #   #    # 
#    #####  ### #     #       #####  ######   #   #    # 




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






### SIRV

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
        Low  = G[hl].vs["S2V_l"][0]
        High = G[hl].vs["S2V_h"][0]
        delta = Low - High
        G[hl].vs["chances_to_vaccinate"] = 1 - Low + (1-PB)*delta    # if this is static, I calculate it only once, and then I use it in the update function


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
            Low  = g_h.vs["S2V_l"][0]
            High = g_h.vs["S2V_h"][0]
            delta = Low - High
            g_b.vs["chances_to_vaccinate"] = 1 - Low + (1-PB)*delta

        S2I = g_h.vs["S2I"][0] # it's the same for everyone, I put it here to save time, faster than accessing the dictionary every time
        I2R = g_h.vs["I2R"][0] # it's the same for everyone, I put it here to save time, faster than accessing the dictionary every time
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(len(g_h.vs)))

        # first, people who wants to vaccinate will try to do so.
        for i,vertex in enumerate(g_h.vs):
            # only a susceptible person can vaccinate
            if vertex["health_status"]==1:
                if(rr[i]<g_h.vs[i]["chances_to_vaccinate"]):
                    vertex['health_status']=4
                else:
                    vertex['health_status']=1
        # these are people who are trying to be vaccinated. if more than VA * N people are trying to be vaccinated, only VA * N will be vaccinated, they are picked at random

        indx_vacc = np.where(np.array(g_h.vs["health_status"])==4)[0]
        if(len(indx_vacc)>g_h.vs["VA"][0]*len(g_h.vs)):
            # if more than VA * N people are trying to be vaccinated, only VA * N will be vaccinated, they are picked at random
            #select a random sample from indx_vacc of len = len(indx_vacc) - VA * N # they will fail to get vaccinated
                                                                                    # and their status will be set back to 1
            L = len(indx_vacc) - int(g_h.vs["VA"][0]*len(g_h.vs))
            for gg in np.random.choice(indx_vacc, L, replace = False):
                g_h.vs[gg]['health_status']=1
        # now, a normal update of the SIR model
        for i,vertex in enumerate(g_h.vs):
            if vertex["health_status"]==1:
                h_n = np.sum(np.array(g_h.vs[g_h.neighbors(i)]["health_status"])==2)

                if(rr[i]<1-np.power(1-S2I,h_n)):
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
            elif vertex["health_status"]==4:
                vertex['temp_status']=4
            else:
                input('ERROR! This node had an health state that was unexpected!!! in evolve_SIRB_beta')

        # update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"]=vertex['temp_status']


# dictionary with the list of rules and the corresponding function to initialize them. to use it call init_single_rule[P_dyn_i["func"]](P_dyn_i, G)

update_dictionary = {
    "SIRB_beta" : update_SIRB_beta,
    "SIRV" : update_SIRV
}

init_single_rule = {
    "SIRB_beta": init_SIRB_beta,
    "SIRV" : init_SIRV
}


# Ask Arne how to make the code expandable without expanding the code. I want users to be able to add new rules without having to modify the code.
# Update schedule?
    # how to do the savings in case of evolutionary game theory? (several waves of infection)

# Same models with different communication regimes (121, 12M, M21)




# Ability to read csv files of graphs
# Temporal network 
# Ability to add new save functions



