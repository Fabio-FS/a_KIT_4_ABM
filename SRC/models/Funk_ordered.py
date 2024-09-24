import sys
sys.path.append('..')
from utilities.IC import *

def calc_TA(graph,p):
    susc_nodes = np.sum(np.array(graph.vs["health_status"])==1)
    if susc_nodes == 0:
        TotAwa = np.nan
        beta_susc = np.nan
    else:
        TotAwa = 0
        beta_susc = 0
        for i,vertex in enumerate(graph.vs):
            if vertex["health_status"] == 1:
                TotAwa += (p**vertex["awareness"])/susc_nodes
                beta_susc += vertex["beta"]/susc_nodes
    return TotAwa, beta_susc

def dbb(vertex,rate,rnd_number):
    if vertex["health_status"] == 2 and rnd_number < rate:
        return 0
    return vertex["next_awareness"]

def social_influence_Funk(vertex,v_index,g_b,rate,random_numbers):
    nghbr_awareness = np.array(g_b.vs[g_b.neighbors(v_index)]["awareness"])   #neighbor's current awareness
    nghbr_awareness = np.sort(nghbr_awareness[nghbr_awareness < vertex["next_awareness"] - 1])

    for nb in nghbr_awareness:
        rnd_number, random_numbers = random_numbers[:1], random_numbers[1:]
        if rnd_number < rate:
            return nb+1, random_numbers #the neighbor awarenesses are sorted. So if one is accepted, the other's have no chance to spread (bc they're worse)
        
    return vertex["next_awareness"], random_numbers
            

def decay(current, rate, rnd_number, max_value):
        print("this is current", current)
        print(current,max_value,rnd_number,rate)
        if current < max_value and rnd_number < rate:
            return current+1
        return current

def update_Funk_ordered(G,rule,global_var):
    """this simulates an SIR model, almost the complete model from Funk 2009.

    Awareness could be handled implicitly, but is carried with explicitly for precision reasons


    Awareness is implicitly taken care and is used to calculate agent's personal betas.
    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer] and in rule
    """

    g_h = G[rule["hl"]]
    g_b = G[rule["bl"]]

    lamda = rule["lambda"]
    beta0 = rule["beta0"]
    p = rule["p"]
    alpha = rule["alpha"]
    omega = rule["omega"]
    max_awareness = rule["max_awareness"]
    max_behavior = rule["max_behavior"]

    rr_1 = np.random.uniform(low=0, high=1, size=g_b.vcount()) < omega  #source
    rr_2 = np.random.uniform(low=0, high=1, size=g_b.vcount()) < lamda   #decay
    rr_3 = np.random.uniform(low = 0, high=1, size = g_b.ecount() *(2 - g_b.is_directed())    )     #awareness spread, random number for each edge

    #calculate the total awareness defined in Funk 2009
    TotAwa, PB_susc = calc_TA(g_h,p)
    G[0].vs["TotAwa"] = TotAwa
    G[0].vs["PB_susc"] = PB_susc

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    #source: infected agents can become aware   (Disease-Behavior-Bridge)
    g_b.vs["next_awareness"] *= 1- (          rr_1             *                 (np.array(g_h.vs["health_status"]) == 2)                 )
    #                                       array 1                                         array 2
    #                               This is zero, where both individual arrays are True
    #                               rr_1 is True where a random value was smaller than omega
    #                               the second array is True where a node is infected
    #                                so this is "if node is infected and random < omega: multiply next_awareness with 0 else 1"

    #spread: aware agents spread awareness       (behavior model)
    for i,vertex in enumerate(g_b.vs):
        #all agents who are more aware than me, will try to spread their awareness to me
        vertex["next_awareness"],rr_3 = social_influence_Funk(vertex,i,g_b,alpha,rr_3)
    
    #decay: aware agents lose awareness          (behavior model)
    g_b.vs["next_awareness"] += (rr_2)*(np.array(g_b.vs["next_awareness"]) < max_awareness)

    #calculate new beta
    g_b.vs["next_beta"] = (1-max_behavior)*beta0 + max_behavior*(1- p ** np.array(g_b.vs["next_awareness"])) * beta0

    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas and awarenesses
    for i,vertex in enumerate(g_b.vs):
        vertex["awareness"] = vertex["next_awareness"]
        vertex["beta"]      = vertex["next_beta"]

    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate update of health status

    # check if there are infected nodes, if not, skip the update
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
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
                input('ERROR! This node had a health state that was unexpected!!! in update_Funk_reworked')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        #f_infected = np.mean(np.array(G[0].vs["health_status"])==2)
        #G[0].vs["I_peak"] = max(G[0].vs[0]["I_peak"], f_infected)
        #print(G[0].vs["I_peak"])
        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

def init_Funk_ordered(P_dyn, G,global_var):

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    G[hl].vs["I2R"]   = P_dyn["I2R"]                                                            # for each node sets the gamma
    G[bl].vs["PB_susc"] = P_dyn["HEALTH"]["beta0"]

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    #N_infected = np.sum(np.array(G[hl].vs["health_status"])==2)
    #print(N_infected)

    #set_continuous_initial_condition(P_dyn["BEHAVIOR"]["IC"], "personal_beta", G[bl])
    #set_continuous_iniitial_condition is not needed. Every agent starts with the same personal beta of beta0
    #I could include beta0 in the IC flag, but then I would have to write it twice. So i just do it here
    G[bl].vs["beta"]      = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])
    G[bl].vs["next_beta"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])

    max_awareness = P_dyn["BEHAVIOR"]["max_awareness"]

    G[bl].vs["awareness"] = np.full( shape=len(G[bl].vs), fill_value = max_awareness)    
    G[bl].vs["next_awareness"] = np.full( shape=len(G[bl].vs), fill_value = max_awareness)
    #initializing awareness as a high value (e.g. 99)
    #awareness goes from 0 to max_awareness. beta is calculated differently than in paper.
    
    p = P_dyn["BEHAVIOR"]["p"]
    max_behavior = P_dyn["BEHAVIOR"].get("max_behavior")   #ranging from 0 to 1. If 1, agents can completely reduce their beta to 0. No chance of getting infected.
    #If 0.5, the best they can do is: beta_i = 0.5 * beta0
    if max_behavior == None:
        max_behavior = 1

    G[bl].vs["TotAwa"] = p**max_awareness  
    G[hl].vs["I_peak"] = 0   

    rule  = {
        'rule': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'lambda' : P_dyn["BEHAVIOR"]["lambda"],
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'p' : p,
        'alpha' : P_dyn["BEHAVIOR"]["alpha"],
        'omega' : P_dyn["BEHAVIOR"]["omega"],
        'max_awareness' : max_awareness,
        "max_behavior" : max_behavior
    }
    return rule

def init_model(update_fct_dict, init_fct_dict):
    update_fct_dict["Funk_ordered"] = update_Funk_ordered
    init_fct_dict["Funk_ordered"] = init_Funk_ordered










