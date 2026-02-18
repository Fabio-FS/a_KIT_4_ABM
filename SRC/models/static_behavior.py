import sys
sys.path.append('..')
from utilities.IC import *

def update_static(G,rule, global_var):

    """this simulates an SIR model with static behavior.
    behavior as a fixed chance to protect in each timestep.

    Each agent can either protect or not.
    They base their decision on:    1. own previous decision
                                    2. neighbor's previous decision
                                    3. global infected

    Agents are influenced by the share of adopting neighbors

    The parameter mu controls how much agents weigh their own behavior versus the social influence.
    If mu is one agents are not influenced by others, if mu is zero agents don't stick to their opinion.

    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer] and in rule
    """

    g_h = G[rule["hl"]]
    g_b = G[rule["bl"]]

    beta0 = rule["beta0"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        global_var.stop_condition = True
        #no return here.
        #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas and awarenesses
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    g_b.vs["behavior"] = (np.array(g_b.vs["probability"]) > rr1).astype(int).tolist()
    g_b.vs["beta"] = ((1-max_behavior*np.array(g_b.vs["behavior"]))*beta0).tolist()

    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate update of health status

    if(N_infected>0):

        I2R = g_h.vs["I2R"][0] # it's the same for everyone, I put it here to save time, faster than accessing the dictionary every time
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(len(g_h.vs)))
        for i,vertex in enumerate(g_h.vs):
            # the updated health status is not overwritten immediately into health_status, but it is stored in next_health to avoid updating the health status of a
            # subsequent node with an updated value. this keeps the update syncronous.
             
            # for each node i, if it is susceptible, check if it gets infected:
            if vertex["health_status"]==1:
                # number of neighbors that are infected. they are already in the right layer, so no need to filter them
                inf_nghbrs = np.sum(np.array(g_h.vs[g_h.neighbors(i)]["health_status"])==2)

                if(rr[i] < 1 - np.power (1 - g_b.vs[i]["beta"],inf_nghbrs)):
                    #      1 -          (1 - beta) ^ I
                    vertex['next_health'] = 2
                #else: vertex['next_health'] = 1     #written here for clarity
            elif(vertex["health_status"] == 2):
                if(rr[i]<I2R):
                    vertex['next_health'] = 3
                #else: vertex['next_health'] = 2      #written here for clarity
            elif vertex["health_status"] == 3:        #condition is needed for the error message below to work
                pass
                #vertex['next_health'] = 3            #written here for clarity
            else:
                input('ERROR! This node had a health state that was unexpected!!! in update_downward')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))


def init_static(P_dyn, G, global_var):

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted
    

    G[hl].vs["I2R"]   = P_dyn["HEALTH"]["I2R"]                                                            # for each node sets the gamma
    G[bl].vs["PB_susc"] = P_dyn["HEALTH"]["beta0"]

    # sets the initial condition for each node
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    G[bl].vs["beta"] = (np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])).astype(float).tolist()
    G[bl].vs["probability"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["BEHAVIOR"]["static_probability"])
    G[hl].vs["I_peak"] = 0
    G[bl].vs["behavior"] = 0

    rule  = {
        'func': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }
    return rule

def init_model(update_fct_dict, init_fct_dict):
    update_fct_dict["static"] = update_static
    init_fct_dict["static"] = init_static