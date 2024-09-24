import sys
sys.path.append('..')
from utilities.IC import *


#to eventually be evaluated in a 2 island setting

def calc_TA(graph,p,idcs=None):
    #through idcs, you can input the indices of vertices that should be evauluated. If None, all vertices will be evaluated

    if idcs is None:
        idcs = range(graph.vcount())

    susc_nodes = np.sum(np.array(graph.vs[idcs]["health_status"])==1)
    if susc_nodes == 0:
        TotAwa = np.nan
        beta_susc = np.nan
    else:
        TotAwa = 0
        beta_susc = 0

        for vertex in graph.vs[idcs]:
            if vertex["health_status"] == 1:
                TotAwa += (p**vertex["awareness"])/susc_nodes
                beta_susc += vertex["beta"]/susc_nodes

    #print(graph.vs["awareness"])
    #print(graph.vs["beta"])
    #print(graph.vs[idcs]["awareness"])
    #print(TotAwa,beta_susc)
    #input()
    return TotAwa, beta_susc

def update_Funk_2i_stoch(G,rule,global_var):
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


    #------------------------------------------------------------------------------------------------------------------------------#
    #first: update awareness
    type_event = np.random.choice([0,1,2], size = g_b.vcount(), p = [alpha/(alpha+omega+lamda),omega/(alpha+omega+lamda),lamda/(alpha+omega+lamda)])

    for i,evt in enumerate(type_event):
        if evt == 0:   #spread
            j = np.random.choice(g_b.neighbors(i))
            g_b.vs[i]["next_awareness"] = min(g_b.vs[i]["next_awareness"],g_b.vs[j]["awareness"]+1)
        elif evt == 1:   #source
            g_b.vs[i]["next_awareness"] = g_b.vs[i]["next_awareness"] * (g_b.vs["health_status"] != 2)
        else:   #decay
            g_b.vs[i]["next_awareness"] += 1

    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas
    g_b.vs["awareness"] = g_b.vs["next_awareness"]
    #calculate the total awareness defined in Funk 2009
    TotAwa, PB_susc = calc_TA(g_h,p)
    G[0].vs["TotAwa"] = TotAwa
    G[0].vs["PB_susc"] = PB_susc
    g_b.vs["beta"] = beta0 * (1 - max_behavior * p ** np.array(g_b.vs["awareness"]))

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
                    vertex['next_health'] = 1
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

def update_Funk_2i_det(G,rule,global_var):
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
    network_is_island = rule["network_is_island"]


    #------------------------------------------------------------------------------------------------------------------------------#
    #first: update awareness
    #type_event = np.random.choice([0,1,2], size = g_b.vcount(), p = [alpha/(alpha+omega+lamda),omega/(alpha+omega+lamda),lamda/(alpha+omega+lamda)])

    for i,vertex in enumerate(g_b.vs):
        j = np.random.choice(g_b.neighbors(i))
        #print(g_b.vs[j]["awareness"], g_b.vs[i]["health_status"])
        g_b.vs[i]["next_awareness"] = (alpha/(alpha+omega+lamda) * min(g_b.vs[i]["next_awareness"],g_b.vs[j]["awareness"]+1) +
                                       omega/(alpha+omega+lamda) * g_b.vs[i]["next_awareness"] * (g_b.vs[i]["health_status"] != 2) +
                                       lamda/(alpha+omega+lamda) * g_b.vs[i]["next_awareness"] + 1)
        #print( g_b.vs[i]["next_awareness"],  g_b.vs[i]["awareness"])
    #input()

    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas
    g_b.vs["awareness"] = g_b.vs["next_awareness"]
    #calculate the total awareness defined in Funk 2009
    TotAwa, PB_susc = calc_TA(g_h,p)
    G[0].vs["TotAwa"] = TotAwa
    G[0].vs["PB_susc"] = PB_susc
    g_b.vs["beta"] = beta0 * (1 - max_behavior * p ** np.array(g_b.vs["awareness"]))

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
                    vertex['next_health'] = 1
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

    if network_is_island:
        g_h.vs["I0"] = np.mean(np.array(g_h.vs[np.where(np.array(g_h.vs["membership"])==0)[0]]["health_status"]) == 2)
        g_h.vs["I1"] = np.mean(np.array(g_h.vs[np.where(np.array(g_h.vs["membership"])==1)[0]]["health_status"]) == 2)
        g_b.vs["TA0"] = calc_TA(g_b, p, idcs = np.where(np.array(g_h.vs["membership"])==0)[0])[0]
        g_b.vs["TA1"] = calc_TA(g_b, p, idcs = np.where(np.array(g_h.vs["membership"])==1)[0])[0]
        #print(g_b.vs["TA0"][0],g_b.vs["TA1"][0])
        #input()

def init_Funk_2i(P_dyn, G,global_var):

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    G[hl].vs["I2R"]   = P_dyn["I2R"]                                                            # for each node sets the gamma
    G[bl].vs["PB_susc"] = P_dyn["HEALTH"]["beta0"]

    # for each node sets the initial condition
    
    if(P_dyn["HEALTH"]["IC"]["type"] == "random"):
        #pick N_pat0 random patient zeros, but pick them only from group 2
        G[hl].vs["health_status"]=1
        # select N_pat_zero nodes at random and set them to 2:
        for gg in rn.sample(range(G[hl].vcount()//2,G[hl].vcount()), P_dyn["HEALTH"]["IC"]["N_pat_zero"]):
            G[hl].vs[gg]["health_status"] = 2
    else:
        set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])

    p = P_dyn["BEHAVIOR"]["p"]
    max_behavior = P_dyn["BEHAVIOR"].get("max_behavior")   #ranging from 0 to 1. If 1, agents can completely reduce their beta to 0. No chance of getting infected.
    #If 0.5, the best they can do is: beta_i = 0.5 * beta0
    #If 1, the best they can do is: beta_i = 0, If 0, they cannot alter beta0 at all
    if max_behavior == None:
        max_behavior = 1

    G[hl].vs["next_health"] = G[hl].vs["health_status"]

    G[bl].vs["beta"]      = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])

    max_awareness = P_dyn["BEHAVIOR"]["max_awareness"]

    G[bl].vs["awareness"] = np.full( shape=len(G[bl].vs), fill_value = max_awareness)
    G[bl].vs["next_awareness"] = np.full( shape=len(G[bl].vs), fill_value = max_awareness)
    #initializing awareness as a high value (e.g. 99)
    #awareness goes from 0 to max_awareness. beta is calculated differently than in paper.

    try:
        G[hl].vs["I0"] = 0
        G[hl].vs["I1"] = np.mean(np.array(G[hl].vs[np.where(np.array(G[hl].vs["membership"])==1)[0]]["health_status"]) == 2)
        G[bl].vs["TA0"] = calc_TA(G[bl], p, idcs = np.where(np.array(G[bl].vs["membership"])==0)[0])[0]
        G[bl].vs["TA1"] = calc_TA(G[bl], p, idcs = np.where(np.array(G[bl].vs["membership"])==1)[0])[0]
        network_is_island = True
    except KeyError:
        network_is_island = False

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
        "max_behavior" : max_behavior,
        "network_is_island" : network_is_island
    }
    return rule

def init_model(update_fct_dict, init_fct_dict):
    update_fct_dict["Funk_2i_stoch"] = update_Funk_2i_stoch
    init_fct_dict["Funk_2i_stoch"] = init_Funk_2i
    update_fct_dict["Funk_2i_det"] = update_Funk_2i_det
    init_fct_dict["Funk_2i_det"] = init_Funk_2i










