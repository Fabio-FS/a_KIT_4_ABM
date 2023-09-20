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
    #print(IC)
    if(IC["type"] == "random"):
        g.vs[name]=1
        # select N_pat_zero nodes at random and set them to 2:
        for gg in rn.sample(range(len(g.vs)), IC["N_pat_zero"]):
            g.vs[gg][name] = 2
    elif(IC["type"] == "central_geometric"):
        g.vs[name]=1
        # select N_pat_zero nodes at random and set them to 2:
        gg = np.floor(len(g.vs)/2)
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
                                                                                "homophily" : { "Flag" : "True",
                                                                                                "target" : -1,
                                                                                                "steps" : 10000,
                                                                                                "precision" : 0.01}}                                                                            }
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
        delta = Low - High
        G[hl].vs["chances_to_stay_healthy"] = 1 - Low + (1-PB)*delta    # if this is static, I calculate it only once, and then I use it in the update function


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
            Low  = g_h.vs["S2I_l"][0]
            High = g_h.vs["S2I_h"][0]
            delta = Low - High
            g_b.vs["chances_to_stay_healthy"] = 1 - Low + (1-PB)*delta


        I2R = g_h.vs["I2R"][0] # it's the same for everyone, I put it here to save time, faster than accessing the dictionary every time
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(len(g_h.vs)))
        for i,vertex in enumerate(g_h.vs):
            # the updated health status is not overwritten immediately into health_status, but it is stored in temp_status to avoid updating the health status of a
            # subsequent node with an updated value. this keeps the update syncronous.
             
            # for each node i, if it is susceptible, check if it gets infected:
            if vertex["health_status"]==1:
                # number of neighbors that are infected. they are already in the right layer, so no need to filter them
                
                #print(np.array(g.vs[g.neighbors(i)]["health_status"])==2)

                h_n = np.sum(np.array(g_h.vs[g_h.neighbors(i)]["health_status"])==2)

                if(rr[i]<np.power(g_h.vs[i]["chances_to_stay_healthy"],h_n)):
                    vertex['temp_status']=1
                else:
                    vertex['temp_status']=2
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
            vertex["health_status"]=vertex['temp_status']

#    #####  ### ######                                           
#   #     #  #  #     #        ####    ##   #    # #    #   ##   
#   #        #  #     #       #    #  #  #  ##  ## ##  ##  #  #  
#    #####   #  ######  ##### #      #    # # ## # # ## # #    # 
#         #  #  #   #         #  ### ###### #    # #    # ###### 
#   #     #  #  #    #        #    # #    # #    # #    # #    # 
#    #####  ### #     #        ####  #    # #    # #    # #    # 
#                                                                


def init_SIRB_gamma(P_dyn, G):

    hl = P_dyn["home_layer"]                # home layer where the dynamic is imprinted
    rule = P_dyn["func"]                    # name of the dynamic

    # imprint on graph:

    G[hl].vs["S2I"] = P_dyn["S2I"]                                          # for each node sets beta
    G[hl].vs["I2R_u"] = P_dyn["I2R_u"]                                      # for each node sets the higher gamma
    G[hl].vs["I2R_l"] = P_dyn["I2R_l"]                                      # for each node sets the lower gamma
    G[hl].vs["read_SIRB_gamma_behavior_from"] = P_dyn["layer_behavior"]      # for each node sets the layer from which to read the behavior
    G[hl].vs["read_as_beh_for_SIRB_gamma"] = P_dyn["name_behavior"]          # for each node sets the name of the behavior
    G[hl].vs["name_for_health_status"] = P_dyn["name"]
    # for each node sets the initial condition

    set_disease_initial_condition(P_dyn["IC"], "health_status", G[hl])
    return hl, rule                                                                                                                    
def update_SIRB_gamma(G, layer):

    """this simulates the classical SIR model, where the probability of being infected depends ONLY on the behavior of the susceptible individuals
    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer]
    """

    # check if there are infected nodes, if not, skip the update
    g = G[layer]
    name = g.vs["name_for_health_status"][0]
    List_index_infected = np.where(np.array(g.vs[name])==2)[0]
    
    if(len(List_index_infected))>0:

        one_minus_beta = 1 - G[layer].vs[0]["S2I"]                  # it's the same for everyone, precalculated to save time

        # for speed: I precalculate the personal gamma ONLY of those who are infected.
        gamma_low = G[layer].vs["I2R_l"][0]
        gamma_high = G[layer].vs["I2R_u"][0]


        read_SIRB_gamma_behavior_from = G[layer].vs["read_SIRB_gamma_behavior_from"][0]
        read_as_beh_for_SIRB_gamma = G[layer].vs["read_as_beh_for_SIRB_gamma"][0]

        # to improve speed further calculate this only for the infected nodes
        PB = np.array(G[read_SIRB_gamma_behavior_from].vs[read_as_beh_for_SIRB_gamma])
        Pgamma = gamma_low + (1-PB)*(gamma_high-gamma_low)

        # generate a random number for each node
        # to improve speed further generate this only for Susceptibles and infected nodes 
        rr=np.random.uniform(low=0, high=1, size=(len(g.vs)))
        for i,vertex in enumerate(g.vs):
            # the updated health status is not overwritten immediately into health_status, but it is stored in temp_status to avoid updating the health status of a
            # subsequent node with an updated value. this keeps the update syncronous.
             
            # for each node i, if it is susceptible, check if it gets infected:
            if vertex[name]==1:
                # number of neighbors that are infected. they are already in the right layer, so no need to filter them
                h_n = np.sum(np.array(g.vs[g.neighbors(i)][name])==2)

                if(rr[i]<np.power(one_minus_beta,h_n)):
                    vertex['temp_status']=1
                else:
                    vertex['temp_status']=2
            elif(vertex[name]==2):
                if(rr[i]<Pgamma[i]):
                    vertex['temp_status']=3
                else:
                    vertex['temp_status']=2
            elif vertex[name]==3:
                vertex['temp_status']=3
            else:
                input('ERROR! This node had an health state that was unexpected!!! in evolve_SIRB_gamma')

        # update the health status of the nodes
        for vertex in g.vs:
            vertex[name]=vertex['temp_status']



#    #####  ### ######                            
#   #     #  #  #     #       #####   ##   #    # 
#   #        #  #     #         #    #  #  #    # 
#    #####   #  ######  #####   #   #    # #    # 
#         #  #  #   #           #   ###### #    # 
#   #     #  #  #    #          #   #    # #    # 
#    #####  ### #     #         #   #    #  ####  
#                                                 

def init_SIRB_tau(P_dyn, G):

    hl = P_dyn["home_layer"]                # home layer where the dynamic is imprinted
    rule = P_dyn["func"]                    # name of the dynamic

    # imprint on graph:


    # TAU is the inverse of gamma. 

    G[hl].vs["S2I"] = P_dyn["S2I"]                                          # for each node sets beta
    G[hl].vs["tau_u"] = P_dyn["tau_u"]                                      # for each node sets the higher gamma
    G[hl].vs["tau_l"] = P_dyn["tau_l"]                                      # for each node sets the lower gamma
    G[hl].vs["read_SIRB_gamma_behavior_from"] = P_dyn["layer_behavior"]      # for each node sets the layer from which to read the behavior
    G[hl].vs["read_as_beh_for_SIRB_gamma"] = P_dyn["name_behavior"]          # for each node sets the name of the behavior
    G[hl].vs["name_for_health_status"] = P_dyn["name"]
    # for each node sets the initial condition

    set_disease_initial_condition(P_dyn["IC"], P_dyn["name"], G[hl])
    return hl, rule
                                                                                                                    
def update_SIRB_tau(G, layer):

    """this simulates the classical SIR model, where the probability of being infected depends ONLY on the behavior of the susceptible individuals
    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer]
    """

    # check if there are infected nodes, if not, skip the update
    g = G[layer]
    name = g.vs["name_for_health_status"][0]
    List_index_infected = np.where(np.array(g.vs[name])==2)[0]
    
    if(len(List_index_infected))>0:

        one_minus_beta = 1 - G[layer].vs[0]["S2I"]                  # it's the same for everyone, precalculated to save time
        

        # for speed: I precalculate the personal gamma ONLY of those who are infected.
        tau_low = G[layer].vs["tau_l"][0]
        tau_high = G[layer].vs["tau_u"][0]


        # to improve speed further calculate this only for the infected nodes
        PB = np.array(G[G[layer].vs["read_SIRB_tau_behavior_from"]].vs["read_as_beh_for_SIRB_tau"])
        Pgamma = 1/ ( tau_high - (1-PB)*(tau_high-tau_low)) # high tau is bad, low tau is good

        # generate a random number for each node
        # to improve speed further generate this only for Susceptibles and infected nodes 
        rr=np.random.uniform(low=0, high=1, size=(len(g.vs)))
        for i,vertex in enumerate(g.vs):
            # the updated health status is not overwritten immediately into health_status, but it is stored in temp_status to avoid updating the health status of a
            # subsequent node with an updated value. this keeps the update syncronous.
             
            # for each node i, if it is susceptible, check if it gets infected:
            if vertex[name]==1:
                # number of neighbors that are infected. they are already in the right layer, so no need to filter them
                h_n = np.sum(np.array(g.vs[g.neighbors(i)][name])==2)

                if(rr[i]<np.power(one_minus_beta,h_n)):
                    vertex['temp_status']=1
                else:
                    vertex['temp_status']=2
            elif(vertex[name]==2):
                if(rr[i]<Pgamma[i]):
                    
                    vertex['temp_status']=3
                else:
                    
                    vertex['temp_status']=2
            elif vertex[name]==3:
                vertex['temp_status']=3
            else:
                input('ERROR! This node had an health state that was unexpected!!! in evolve_SIRB_tau')

        # update the health status of the nodes
        for vertex in g.vs:
            vertex[name]=vertex['temp_status']





#    #####                                    ######                                              
#   #     # #####   ##   ##### #  ####        #     # ###### #    #   ##   #    # #  ####  #####  
#   #         #    #  #    #   # #    #       #     # #      #    #  #  #  #    # # #    # #    # 
#    #####    #   #    #   #   # #      ##### ######  #####  ###### #    # #    # # #    # #    # 
#         #   #   ######   #   # #            #     # #      #    # ###### #    # # #    # #####  
#   #     #   #   #    #   #   # #    #       #     # #      #    # #    #  #  #  # #    # #   #  
#    #####    #   #    #   #   #  ####        ######  ###### #    # #    #   ##   #  ####  #    # 
#                                                                                                 


def init_Static_Behavior(P_dyn, G):                     # This function is used to initialize a variable (usually behavior) that would not be updated during the simulation

    hl = P_dyn["home_layer"]                            # home layer where the dynamic is imprinted

    #print(G[hl])
    # for each node sets the initial condition
    set_continuous_initial_condition(P_dyn["IC"], P_dyn["name"], G[hl])

    return None, None                                   # Since this is not a dynamic, it returns None, None
def update_Static_Behavior(G, layer):
    print("if you read this there is a bug. the code shouldn't have called this function, since there is no dynamic associated to this rule")
    pass

                        


def init_Average(P_dyn, G):

    hl = P_dyn["home_layer"]                # home layer where the dynamic is imprinted
    rule = P_dyn["func"]                    # name of the dynamic

    # imprint on graph:

    G[hl].vs["mu"] = P_dyn["mu"]                       # for each node sets the mu
    
    G[hl].vs["name_for_behavior_status"] = P_dyn["name"]
    # for each node sets the initial condition

    set_continuous_initial_condition(P_dyn["IC"], P_dyn["name"], G[hl])

    return hl, rule
def update_Average(G, layer):
    # for each node, calculate the average of the neighbors and update the value of the node
    g = G[layer]
    name = g.vs["name_for_behavior_status"][0]
    for vertex in g.vs:

        list_of_neighbors = vertex.neighbors()

        # select the value of "behavior_status" of a random neighbor
        if list_of_neighbors:
            random_value = g.vs[np.random.choice(list_of_neighbors)][name]
        else:
            random_value = np.random.random()

        # update the value of the node new = old*(1-mu) + random_value*mu

        vertex["temp_status"] = vertex[name]*(1-vertex["mu"]) +  random_value*vertex["mu"]

    # synchronize the update
    for vertex in g.vs:
        vertex[name] = vertex["temp_status"]





# dictionary with the list of rules and the corresponding function to initialize them. to use it call init_single_rule[P_dyn_i["func"]](P_dyn_i, G)



update_dictionary = {
    "SIRB_beta" : update_SIRB_beta,
    "SIRB_gamma" : update_SIRB_gamma,
    "SIRB_tau" : update_SIRB_tau,
    "Average" : update_Average,
    "Static_Behavior" : update_Static_Behavior
}

init_single_rule = {
    "SIRB_beta": init_SIRB_beta,
    "SIRB_gamma": init_SIRB_gamma,
    "SIRB_tau": init_SIRB_tau,
    "Average" : init_Average,
    "Static_Behavior" : init_Static_Behavior
}
