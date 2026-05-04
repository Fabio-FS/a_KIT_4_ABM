import sys
sys.path.append('..')
from utilities.IC import *
from utilities.model_helpers import *

def update_static_probability(G,rule, global_var):

    """this simulates an SIR model with static behavior, i.e. the probability to protect doesn't change.
    behavior as a fixed chance to protect in each timestep.

    Each agent can either protect or not.

    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer] and in rule
    """

    g_b = G[rule["bl"]]

    beta0 = rule["beta0"]
    max_behavior = rule["max_behavior"]

    health_status = global_var.health_status

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum( health_status == 2)
    if N_infected == 0:
        N_protecting = np.sum(global_var.behavior == 1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    global_var.behavior = (global_var.probability > rr1).astype(int)
    global_var.beta = ((1-max_behavior*global_var.behavior)*beta0)

    #------------------------------------------------------------------------------------------------------------------------------#
    #third: update the health status of the nodes
    global_var.health_status = update_health_SIR(N_infected, global_var)

    global_var.I_peak = max(global_var.I_peak  ,  float(np.mean(health_status==2)))

def update_static(G,rule, global_var):

    """this simulates an SIR model with static behavior, i.e. the behavior doesn't change.

    Each agent is either protecting or not.

    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer] and in rule
    """

    health_status = global_var.health_status

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum( health_status == 2)
    if N_infected == 0:
        global_var.stop_condition = True
        #no return here.
        #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the health status of the nodes
    global_var.health_status = update_health_SIR(N_infected, global_var)

    global_var.I_peak = max(global_var.I_peak  ,  float(np.mean(health_status==2)))

def init_static_probability(P_dyn, G, global_var):

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted
    

    global_var.I2R   = P_dyn["HEALTH"]["I2R"]                                                            # for each node sets the gamma
    #G[bl].vs["PB_susc"] = P_dyn["HEALTH"]["beta0"]

    global_var.h_neighbors = [np.array(G[hl].neighbors(i)) for i in range(G[hl].vcount())]
    global_var.N_nodes_H = G[hl].vcount()

    # sets the initial condition for each node
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])

    global_var.health_status = np.array(G[hl].vs["health_status"])
    global_var.beta = np.full( shape = G[bl].vcount(), fill_value = P_dyn["HEALTH"]["beta0"])
    global_var.probability = np.full( shape = G[bl].vcount(), fill_value = P_dyn["BEHAVIOR"]["static_probability"])
    global_var.I_peak = 0
    global_var.behavior = np.zeros(shape = G[bl].vcount())

    rule  = {
        'func': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }
    return rule

def init_static(P_dyn, G, global_var):

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted
    
    global_var.I2R   = P_dyn["HEALTH"]["I2R"]                                                            # for each node sets the gamma

    global_var.h_neighbors = [np.array(G[hl].neighbors(i)) for i in range(G[hl].vcount())]
    global_var.N_nodes_H = G[hl].vcount()

    # sets the initial condition for each node
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])

    global_var.health_status = np.array(G[hl].vs["health_status"])
    global_var.behavior = np.random.random(size = G[bl].vcount())  <  P_dyn["BEHAVIOR"]["static_probability"]

    if P_dyn["BEHAVIOR"]["IC"]["homophily"]["Flag"] == True and P_dyn["BEHAVIOR"]["static_probability"] != 1 and P_dyn["BEHAVIOR"]["static_probability"] != 0:  #not doing calculations in homogeneous populations
        set_initial_condition(P_dyn["BEHAVIOR"]["IC"],
                              attribute_vector = global_var.behavior,
                              g = G[bl],
                              global_var = global_var)


    global_var.beta = np.where(global_var.behavior == 1,
                               (1 -  P_dyn["BEHAVIOR"]["max_behavior"]) * P_dyn["HEALTH"]["beta0"] ,
                               P_dyn["HEALTH"]["beta0"])
    global_var.I_peak = 0

    rule  = {
        'func': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }
    return rule

def init_model(update_fct_dict, init_fct_dict):
    update_fct_dict["static_prob"] = update_static_probability
    init_fct_dict["static_prob"] = init_static_probability

    update_fct_dict["static"] = update_static
    init_fct_dict["static"] = init_static