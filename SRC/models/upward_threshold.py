import sys
sys.path.append('..')
from utilities.IC import *

def update_upward(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    a_corr = rule["a_corr"]
    pn_thr = rule["pn_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)      
        vertex["probability"] = ( (1 / (1 + np.exp(   - a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                      - a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                      - a_Ni*           (N_infected/g_h.vcount() - Ni_thr))) )  *
                                  (2 / (1 + np.exp(-a_corr *           (N_infected/g_h.vcount()) )) - 1) )
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
                input('ERROR! This node had a health state that was unexpected!!! in update_upward')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

def update_upward_nocorr(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    pn_thr = rule["pn_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)      
        vertex["probability"] = 1 / (1 + np.exp(   - a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                   - a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                   - a_Ni*           (N_infected/g_h.vcount() - Ni_thr)))
    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas and awarenesses
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    g_b.vs["behavior"] = (np.array(g_b.vs["probability"]) > rr1).astype(int).tolist()
    if N_infected == 0:
        g_b.vs["behavior"] = 0 #not completely true to the model but doesn't change any macro prediction and speeds up simulations dramatically (probably)
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
                input('ERROR! This node had a health state that was unexpected!!! in update_upward')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

def update_upw_mov(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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

    #I needed this extra function for making movies
    #the other functions update behavior, update disease, then take snapshot
    #but in movies, I want to show the behavior that will have an impact on the next step
    #therefore I update disease, then behavior.
    #But the first update should still be behavior.
    #So I am doing    behavior  -  snapshot -  disease,behavior  - snapshot - ....

    g_h = G[rule["hl"]]
    g_b = G[rule["bl"]]

    beta0 = rule["beta0"]
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    a_corr = rule["a_corr"]
    pn_thr = rule["pn_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    if global_var.first_tick == False:
        #health update will not be carried out in the very first step
        #that is why I need a special function for movies upw_MOV UPW_MOV
        #so I can record the initial condition
        #------------------------------------------------------------------------------------------------------------------------------#
        #first: calculate update of health status

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
                    input('ERROR! This node had a health state that was unexpected!!! in update_upward')

            #------------------------------------------------------------------------------------------------------------------------------#
            #second: update the health status of the nodes
            for vertex in g_h.vs:
                vertex["health_status"] = vertex['next_health']

            G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    global_var.first_tick = False
    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)
        vertex["probability"] = ( (1 / (1 + np.exp(   - a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                      - a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                      - a_Ni*           (N_infected/g_h.vcount() - Ni_thr))) )  *
                                  (2 / (1 + np.exp(-a_corr *           (N_infected/g_h.vcount()) )) - 1) )
    #------------------------------------------------------------------------------------------------------------------------------#
    #fourth: update the betas and awarenesses
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    g_b.vs["behavior"] = (np.array(g_b.vs["probability"]) > rr1).astype(int).tolist()
    g_b.vs["beta"] = ((1-max_behavior*np.array(g_b.vs["behavior"]))*beta0).tolist()


def update_upward_Heav(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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
    mu = rule["mu"]
    nifactor = rule["nifactor"]
    thr = rule["thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)      
        vertex["probability"] = np.heaviside(mu*protecting_nghbrs + (1 - mu)*vertex["behavior"] + nifactor*N_infected/g_h.vcount() - thr,1)
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
                input('ERROR! This node had a health state that was unexpected!!! in update_upward')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

def update_downward(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    a_corr = rule["a_corr"]
    pn_thr = rule["pn_thr"]
    Bi_thr = rule["Bi_thr"]
    Ni_thr = rule["Ni_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas
    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)
        vertex["probability"] = ( (1 / (1 + np.exp(     a_pn * mu     * (protecting_nghbrs-pn_thr)
                                                      - a_Bi * (1-mu) * (vertex["behavior"]-Bi_thr) 
                                                      - a_Ni* (N_infected/g_h.vcount()-Ni_thr))) )  *
                                  (2 / (1 + np.exp(-a_corr * (N_infected/g_h.vcount()) )) - 1) )
    
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


def update_downward_nocorr(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    pn_thr = rule["pn_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)      
        vertex["probability"] = 1 / (1 + np.exp(   + a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                   - a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                   - a_Ni*           (N_infected/g_h.vcount() - Ni_thr)))
    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas and awarenesses
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    g_b.vs["behavior"] = (np.array(g_b.vs["probability"]) > rr1).astype(int).tolist()
    if N_infected == 0:
        g_b.vs["behavior"] = 0 #not completely true to the model but doesn't change any macro prediction and speeds up simulations dramatically (probably)
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
                input('ERROR! This node had a health state that was unexpected!!! in update_upward')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

def update_downward_Heav(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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
    mu = rule["mu"]
    nifactor = rule["nifactor"]
    thr = rule["thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)      
        vertex["probability"] = np.heaviside(mu*(1 - protecting_nghbrs) + (1 - mu)*vertex["behavior"] + nifactor*N_infected/g_h.vcount() - thr,1)
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
                input('ERROR! This node had a health state that was unexpected!!! in update_upward')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

def update_dow_mov(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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

    #I needed this extra function for making movies
    #the other functions update behavior, update disease, then take snapshot
    #but in movies, I want to show the behavior that will have an impact on the next step
    #therefore I update disease, then behavior.
    #But the first update should still be behavior.
    #So I am doing    behavior  -  snapshot -  disease,behavior  - snapshot - ....

    g_h = G[rule["hl"]]
    g_b = G[rule["bl"]]

    beta0 = rule["beta0"]
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    a_corr = rule["a_corr"]
    pn_thr = rule["pn_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    if global_var.first_tick == False:
        #------------------------------------------------------------------------------------------------------------------------------#
        #first: calculate update of health status

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
                    input('ERROR! This node had a health state that was unexpected!!! in update_upward')

            #------------------------------------------------------------------------------------------------------------------------------#
            #second: update the health status of the nodes
            for vertex in g_h.vs:
                vertex["health_status"] = vertex['next_health']

            G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    global_var.first_tick = False
    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)
        vertex["probability"] = ( (1 / (1 + np.exp(     a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                      - a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                      - a_Ni*           (N_infected/g_h.vcount() - Ni_thr))) )  *
                                  (2 / (1 + np.exp(-a_corr *           (N_infected/g_h.vcount()) )) - 1) )
    #------------------------------------------------------------------------------------------------------------------------------#
    #fourth: update the betas and awarenesses
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    g_b.vs["behavior"] = (np.array(g_b.vs["probability"]) > rr1).astype(int).tolist()
    g_b.vs["beta"] = ((1-max_behavior*np.array(g_b.vs["behavior"]))*beta0).tolist()

def update_doped(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    pn_thr = rule["pn_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)      
        vertex["probability"] = 1 / (1 + np.exp(     (1 - 2*vertex["herder"]) * a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                   -                            a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                   - a_Ni*           (N_infected/g_h.vcount() - Ni_thr)))
    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas and awarenesses
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    g_b.vs["behavior"] = (np.array(g_b.vs["probability"]) > rr1).astype(int).tolist()
    if N_infected == 0:
        g_b.vs["behavior"] = 0 #not completely true to the model but doesn't change any macro prediction and speeds up simulations dramatically (probably)
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
                input('ERROR! This node had a health state that was unexpected!!! in update_upward')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

def update_mix_3_populations(G, rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    pn_thr = rule["pn_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs[global_var.herder_idcs]):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)      
        vertex["probability"] = 1 / (1 + np.exp(   - a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                   - a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                   - a_Ni*           (N_infected/g_h.vcount() - Ni_thr)))
    for i,vertex in enumerate(g_b.vs[global_var.contrarian_idcs]):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)      
        vertex["probability"] = 1 / (1 + np.exp(   + a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                   - a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                   - a_Ni*           (N_infected/g_h.vcount() - Ni_thr)))
    #no probability update for static agents
    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas and awarenesses
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    g_b.vs["behavior"] = (np.array(g_b.vs["probability"]) > rr1).astype(int).tolist()
    if N_infected == 0:
        g_b.vs["behavior"] = 0 #not completely true to the model but doesn't change any macro prediction and speeds up simulations dramatically (probably)
    g_b.vs["beta"] = ((1-max_behavior*np.array(g_b.vs["behavior"]))*beta0).tolist()

    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate update of health status

    if(N_infected>0):

        I2R = global_var.I2R
        
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
                input('ERROR! This node had a health state that was unexpected!!! in update_upward')

        #------------------------------------------------------------------------------------------------------------------------------#
        #fourth: update the health status of the nodes
        for vertex in g_h.vs:
            vertex["health_status"] = vertex['next_health']

        G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

def update_doped_mov(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to a downward threshold

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

    #I needed this extra function for making movies
    #the other functions update behavior, update disease, then take snapshot
    #but in movies, I want to show the behavior that will have an impact on the next step
    #therefore I update disease, then behavior.
    #But the first update should still be behavior.
    #So I am doing    behavior  -  snapshot -  disease,behavior  - snapshot - ....

    g_h = G[rule["hl"]]
    g_b = G[rule["bl"]]

    beta0 = rule["beta0"]
    mu = rule["mu"]
    a_pn = rule["a_pn"]
    a_Ni = rule["a_Ni"]
    a_Bi = rule["a_Bi"]
    a_corr = rule["a_corr"]
    pn_thr = rule["pn_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(g_b.vs["behavior"])==1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    if global_var.first_tick == False:
        #health update will not be carried out in the very first step
        #that is why I need a special function for movies upw_MOV UPW_MOV
        #so I can record the initial condition
        #------------------------------------------------------------------------------------------------------------------------------#
        #first: calculate update of health status

        if(N_infected>0):

            I2R = g_h.vs["I2R"][0] # it's the same for everyone, I put it here to save time, faster than accessing the dictionary every time
        
            # generate a random number for each node
            rr=np.random.uniform(low=0, high=1, size=(len(g_h.vs)))
            for i,vertex in enumerate(g_h.vs):
                # the updated health status is not overwritten immediately into health_status,
                # the update is stored in next_health, so all health statuses can be updated synchronously
             
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
                    input('ERROR! This node had a health state that was unexpected!!! in update_upward')

            #------------------------------------------------------------------------------------------------------------------------------#
            #second: update the health status of the nodes
            for vertex in g_h.vs:
                vertex["health_status"] = vertex['next_health']

            G[0].vs["I_peak"] = max(G[0].vs["I_peak"][0], np.mean(np.array(G[0].vs["health_status"])==2))

    N_infected = np.sum(np.array(g_h.vs["health_status"])==2)
    global_var.first_tick = False
    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate update of personal betas

    for i,vertex in enumerate(g_b.vs):
        protecting_nghbrs = np.mean(np.array(g_b.vs[g_b.neighbors(i)]["behavior"])==1)
        #vertex["probability"] = vertex["herder"]
        vertex["probability"] = 1 / (1 + np.exp(     (1 - 2*vertex["herder"]) * a_pn * mu     * (protecting_nghbrs       - pn_thr)
                                                   -                            a_Bi * (1-mu) * (vertex["behavior"]      - Bi_thr) 
                                                   - a_Ni*           (N_infected/g_h.vcount() - Ni_thr)))

    #------------------------------------------------------------------------------------------------------------------------------#
    #fourth: update the betas and awarenesses
    rr1 = np.random.uniform(low=0, high=1, size=g_b.vcount())
    g_b.vs["behavior"] = (np.array(g_b.vs["probability"]) > rr1).astype(int).tolist()
    g_b.vs["beta"] = ((1-max_behavior*np.array(g_b.vs["behavior"]))*beta0).tolist()

def init_mix_3_populations(P_dyn, G, global_var):

    global_var.first_tick = True

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    global_var.I2R  = P_dyn["I2R"]                                                       # for each node sets the gamma
    G[bl].vs["PB_susc"] = P_dyn["HEALTH"]["beta0"]

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    G[bl].vs["beta"] = [P_dyn["HEALTH"]["beta0"]]*len(G[bl].vs) #  list(np.full( shape=len(G[bl].vs), fill_value = beta0))
    G[bl].vs["behavior"] = np.full( shape=len(G[bl].vs), fill_value = 0)
    G[bl].vs["next_beta"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])
    G[hl].vs["I_peak"] = 0

    share_herders = P_dyn["BEHAVIOR"]["IC"]["share_herders"]
    share_contrarians = P_dyn["BEHAVIOR"]["IC"]["share_contrarians"]
    share_static = P_dyn["BEHAVIOR"]["IC"]["share_static"]

    if abs(share_contrarians + share_herders + share_static - 1) > 1e-6:
        fixable = False

        if P_dyn["BEHAVIOR"]["IC"]["undefined"] == "static":
            if (share_contrarians + share_herders) < 1:            
                share_static = 1 - (share_contrarians + share_herders)
                fixable = True
        elif P_dyn["BEHAVIOR"]["IC"]["undefined"] == "contrarians":
            if (share_static + share_herders) < 1:            
                share_contrarians = 1 - (share_static + share_herders)
                fixable = True

        if not fixable:
            print("fractions of agents are not adding up to one.")
            exit()

    #assigning all vertices an IRF group.
    remaining_idcs = np.array(range(G[bl].vcount()))
    herder_idcs = np.random.choice(G[bl].vcount(),size=int(share_herders*G[bl].vcount()),replace=False)
    remaining_idcs = np.setdiff1d(remaining_idcs, herder_idcs)
    if remaining_idcs.shape[0] != 0 and share_contrarians != 0:
        contrarian_idcs = np.random.choice(remaining_idcs,size=min(len(remaining_idcs),int(share_contrarians*G[bl].vcount())),replace=False)
        remaining_idcs = np.setdiff1d(remaining_idcs, contrarian_idcs)
    else:
        contrarian_idcs = np.array([])

    #setting up the attribute that the Metropolis algorithm will distribute
    G[bl].vs[herder_idcs]["IRF_group"] = 0
    G[bl].vs[contrarian_idcs]["IRF_group"] = 1
    G[bl].vs[remaining_idcs]["IRF_group"] = 2

    #if homophily flag is True and the population is not homogeneous, REdistribute the attribute "IRF_group"
    if P_dyn["BEHAVIOR"]["IC"]["homophily"]["Flag"] == True and share_herders != 1 and share_contrarians != 1 and share_static != 1:  #not doing calculations in homogeneous populations
        set_initial_condition(P_dyn["BEHAVIOR"]["IC"], 
                                         attribute = "IRF_group",
                                         g = G[bl], 
                                         vector_from_init_fct = np.int32(np.array(G[bl].vs["IRF_group"])))
    
    #after potentially homophilously distriubting the attribute, write the indices to global_var
    global_var.herder_idcs = np.where(np.array(G[bl].vs["IRF_group"]) == 0)[0]
    global_var.contrarian_idcs = np.where(np.array(G[bl].vs["IRF_group"]) == 1)[0]
    global_var.remaining_idcs = np.where(np.array(G[bl].vs["IRF_group"]) == 2)[0]
        
    for vertex in G[bl].vs:
        if vertex["IRF_group"] == 2:
            vertex["probability"] = P_dyn["BEHAVIOR"]["static_probability"]
    #print(G[bl].vs["probability"])
        
    
    #color_dict = {0: (200,0,0), 1: (0,200,0), 2:(0,0,200)}
    #defining the layout
    #mylayout = []
    #for row in range(31):
    #    for col in range(31):
    #        mylayout.append([col,row])
    #G[bl].vs["color"] = [color_dict[i] for i in G[bl].vs["IRF_group"]]
    #from matplotlib import pyplot as plt
    #ig.plot(G[bl])
    #ig.plot(G[bl],target="sample_network.png", layout = mylayout)
    #plt.show()
    #exit()


    rule  = {
        'func': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'mu' : P_dyn["BEHAVIOR"]["mu"],
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'a_pn' : P_dyn["BEHAVIOR"]["a_pn"],
        'a_Ni' : P_dyn["BEHAVIOR"]["a_Ni"],
        'a_Bi' : P_dyn["BEHAVIOR"]["a_Bi"],
        'pn_thr' : P_dyn["BEHAVIOR"]["pn_thr"],
        'Bi_thr' : P_dyn["BEHAVIOR"]["Bi_thr"],
        'Ni_thr' : P_dyn["BEHAVIOR"]["Ni_thr"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }
    return rule



def init_up_down(P_dyn, G,global_var):

    global_var.first_tick = True

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    G[hl].vs["I2R"]   = P_dyn["I2R"]                                                            # for each node sets the gamma
    G[bl].vs["PB_susc"] = P_dyn["HEALTH"]["beta0"]

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    G[bl].vs["beta"] = [P_dyn["HEALTH"]["beta0"]]*len(G[bl].vs) #  list(np.full( shape=len(G[bl].vs), fill_value = beta0))
    G[bl].vs["behavior"] = np.full( shape=len(G[bl].vs), fill_value = 0).tolist()
    G[bl].vs["next_beta"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])
    G[hl].vs["I_peak"] = 0

    a_corr = P_dyn["BEHAVIOR"].get("a_corr",None)

    if P_dyn["func"] == "doped+-" or P_dyn["func"] == "doped_MOV":
        #doped as in semiconductors. Herders and contrarians in one network
        share_herders = P_dyn["BEHAVIOR"]["IC"]["share_herders"]

        herder_idcs = np.random.choice(G[bl].vcount(),size=int(share_herders*G[bl].vcount()),replace=False)
        G[bl].vs["herder"] = False
        G[bl].vs[herder_idcs]["herder"] = True
        if P_dyn["BEHAVIOR"]["IC"]["homophily"]["Flag"] == True and 1 != share_herders != 0:    #not doing calculations in homogeneous populations
            hom = P_dyn["BEHAVIOR"]["IC"]["homophily"]["hom_target"]
            while not (3*share_herders >= 4* hom -0.4 and 
                       3*share_herders >= -4*hom -0.4 and 
                       3*share_herders <= 4*hom + 3.4 and 
                       3*share_herders <= -4*hom + 3.4):
                P_dyn["BEHAVIOR"]["IC"]["homophily"]["hom_target"] = 2* np.random.random() - 1
                hom = P_dyn["BEHAVIOR"]["IC"]["homophily"]["hom_target"]
                #WARNING WARNING WARNING WARNING.
                #This condition is based on a heuristic observation that generally the achievable homophilies satisfy this condition
                #BUT: this is only valid for the network I am currently using, kRRG k=10


            set_initial_condition(P_dyn["BEHAVIOR"]["IC"], 
                                             attribute = "herder",
                                             g = G[bl], 
                                             vector_from_init_fct = np.int32(np.array(G[bl].vs["herder"])))

            G[bl].vs["herder"] = [bool(x == 1) for x in G[bl].vs["herder"]]

    rule  = {
        'func': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'mu' : P_dyn["BEHAVIOR"]["mu"],
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'a_pn' : P_dyn["BEHAVIOR"]["a_pn"],
        'a_Ni' : P_dyn["BEHAVIOR"]["a_Ni"],
        'a_Bi' : P_dyn["BEHAVIOR"]["a_Bi"],
        'a_corr' : a_corr,
        'pn_thr' : P_dyn["BEHAVIOR"]["pn_thr"],
        'Bi_thr' : P_dyn["BEHAVIOR"]["Bi_thr"],
        'Ni_thr' : P_dyn["BEHAVIOR"]["Ni_thr"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }
    return rule

def init_Heav(P_dyn, G,global_var):

    global_var.first_tick = True

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    G[hl].vs["I2R"]   = P_dyn["I2R"]                                                            # for each node sets the gamma
    G[bl].vs["PB_susc"] = P_dyn["HEALTH"]["beta0"]

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    G[bl].vs["beta"] = [P_dyn["HEALTH"]["beta0"]]*len(G[bl].vs) #  list(np.full( shape=len(G[bl].vs), fill_value = beta0))
    G[bl].vs["behavior"] = np.full( shape=len(G[bl].vs), fill_value = 0)
    G[bl].vs["next_beta"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])
    G[hl].vs["I_peak"] = 0

    rule  = {
        'func': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'mu' : P_dyn["BEHAVIOR"]["mu"],
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'nifactor' : P_dyn["BEHAVIOR"]["nifactor"],
        'thr' : P_dyn["BEHAVIOR"]["thr"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }
    return rule

def init_static(P_dyn, G, global_var):

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted
    

    G[hl].vs["I2R"]   = P_dyn["I2R"]                                                            # for each node sets the gamma
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
    update_fct_dict["UPW_CORR"] = update_upward
    update_fct_dict["UPW_MOV"] = update_upw_mov
    update_fct_dict["UPW"] = update_upward_nocorr
    update_fct_dict["UPW_Heav"] = update_upward_Heav  
    update_fct_dict["DOW_CORR"] = update_downward
    update_fct_dict["DOW_MOV"] = update_dow_mov
    update_fct_dict["DOW"] = update_downward_nocorr
    update_fct_dict["DOW_Heav"] = update_downward_Heav
    update_fct_dict["static"] = update_static
    update_fct_dict["doped+-"] = update_doped
    update_fct_dict["doped_MOV"] = update_doped_mov
    update_fct_dict["mix_3"] = update_mix_3_populations
    init_fct_dict["UPW_CORR"] = init_up_down
    init_fct_dict["DOW_CORR"] = init_up_down
    init_fct_dict["static"] = init_static
    init_fct_dict["UPW_MOV"] = init_up_down
    init_fct_dict["DOW_MOV"] = init_up_down
    init_fct_dict["UPW"] = init_up_down
    init_fct_dict["DOW"] = init_up_down
    init_fct_dict["UPW_Heav"] = init_Heav
    init_fct_dict["DOW_Heav"] = init_Heav
    init_fct_dict["doped+-"] = init_up_down
    init_fct_dict["doped_MOV"] = init_up_down
    init_fct_dict["mix_3"] = init_mix_3_populations









