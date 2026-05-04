import sys
sys.path.append('..')
from utilities.IC import *
from utilities.helpers import *
from utilities.model_helpers import *


# ██  ██   █████    ████       ██     ██████   ██████  
# ██  ██   ██  ██   ██ ██     ████      ██     ██      
# ██  ██   ██  ██   ██  ██   ██  ██     ██     ██      
# ██  ██   █████    ██  ██   ██████     ██     ████    
# ██  ██   ██       ██  ██   ██  ██     ██     ██      
# ██  ██   ██       ██ ██    ██  ██     ██     ██      
#  ████    ██       ████     ██  ██     ██     ██████  


def update_upw_dow(G,rule, global_var):

    """this simulates an SIR model
    the behavior model is either an upwards or downwards sloping IRF
    herding vs contrarianism

    Each agent can either protect or not.
    They base their decision on:    1. own previous decision
                                    2. neighbor's previous decision
                                    3. global infected

    The parameter mu controls how much agents weigh their own behavior versus the social influence.
    If mu is one agents are only influenced by others, if mu is one agents don't stick to their behavior

    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer] and in rule
    """

    g_h = G[rule["hl"]]
    g_b = G[rule["bl"]]

    beta0 = rule["beta0"]
    mu = rule["mu"]
    a_B = rule["a_B"]
    a_Ni = rule["a_Ni"]
    BG_thr = rule["BG_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    health_status = global_var.health_status
    behaviors = global_var.behavior

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum( health_status == 2)
    if N_infected == 0:
        N_protecting = np.sum(behaviors == 1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of protection probability

    probability = global_var.functions.calc_protection_probability(behaviors, N_infected, global_var,
                                a_B = a_B, a_Ni = a_Ni, mu = mu,
                                theta = Ni_thr, Bi_thr = Bi_thr, BG_thr = BG_thr)
    
    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas

    global_var.behavior, global_var.beta = update_beta(probability, N_infected, max_behavior, beta0, global_var.N_nodes_H)

    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate and carry out update of health status

    global_var.health_status = update_health_SIR(N_infected, global_var)
    #writing the health status to nodes, because I still use that when saving

    global_var.I_peak = max(global_var.I_peak  ,  float(np.mean(health_status==2)))

def update_upw_dow_health_fixed(G,rule, global_var):

    """this simulates an SIR model
    the behavior model is either an upwards or downwards sloping IRF
    herding vs contrarianism

    Each agent can either protect or not.
    They base their decision on:    1. own previous decision
                                    2. neighbor's previous decision
                                    3. global infected

    The parameter mu controls how much agents weigh their own behavior versus the social influence.
    If mu is one agents are only influenced by others, if mu is one agents don't stick to their behavior

    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer] and in rule
    """

    beta0 = rule["beta0"]
    mu = rule["mu"]
    a_B = rule["a_B"]
    a_Ni = rule["a_Ni"]
    BG_thr = rule["BG_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    health_status = global_var.health_status
    behaviors = global_var.behavior

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum( health_status == 2)
    if N_infected == 0:
        N_protecting = np.sum(behaviors == 1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of protection probability

    probability = global_var.functions.calc_protection_probability(behaviors, N_infected, global_var,
                                a_B = a_B, a_Ni = a_Ni, mu = mu,
                                theta = Ni_thr, Bi_thr = Bi_thr, BG_thr = BG_thr)
    
    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas

    global_var.behavior, global_var.beta = update_beta(probability, N_infected, max_behavior, beta0, global_var.N_nodes_H)

#  ████     ████     ████   
# ██  ██     ██     ██  ██  
# ██         ██     ██      
#  ████      ██      ████   
#     ██     ██         ██  
# ██  ██     ██     ██  ██  
#  ████     ████     ████   

def update_upw_dow_SIS(G,rule, global_var):

    """this simulates an SIS model
    the behavior model is either an upwards or downwards sloping IRF
    herding vs contrarianism

    Each agent can either protect or not.
    They base their decision on:    1. own previous decision
                                    2. neighbor's previous decision
                                    3. global infected

    The parameter mu controls how much agents weigh their own behavior versus the social influence.
    If mu is one agents are only influenced by others, if mu is one agents don't stick to their behavior

    G is the list of graphs-layers,
    layer is the layer where the dynamic is imprinted.
    All the values needed for the simulation are already imprinted in the graph G[layer] and in rule
    """

    g_h = G[rule["hl"]]
    g_b = G[rule["bl"]]

    beta0 = rule["beta0"]
    mu = rule["mu"]
    a_Ni = rule["a_Ni"]
    a_B = rule["a_B"]
    BG_thr = rule["BG_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    health_status = global_var.health_status
    behaviors = global_var.behavior

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum( health_status == 2)
    if N_infected == 0:
        N_protecting = np.sum(behaviors == 1)   #check if there are protecting nodes, if not, everything will be skipped
        if N_protecting == 0:
            global_var.stop_condition = True
            #no return here.
            #batch saving is triggered in the following timestep

    #------------------------------------------------------------------------------------------------------------------------------#
    #first: calculate update of personal betas

    probability = global_var.functions.calc_protection_probability(behaviors, N_infected, global_var,
                                a_B = a_B, a_Ni = a_Ni, mu = mu,
                                theta = Ni_thr, Bi_thr = Bi_thr, BG_thr = BG_thr)
    
    #------------------------------------------------------------------------------------------------------------------------------#
    #second: update the betas

    global_var.behavior, global_var.beta = update_beta(probability, N_infected, max_behavior, beta0, global_var.N_nodes_H)

    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate and carry out update of health status

    global_var.health_status = update_health_SIS(N_infected, global_var)
    g_h.vs["health_status"] = health_status.tolist()
    #writing the health status to nodes, because I still use that when saving

    global_var.I_peak = max(global_var.I_peak  ,  float(np.mean(health_status==2)))


# ██   ██   ████    ██  ██    ████    ██████  
# ███ ███  ██  ██   ██  ██     ██     ██      
# ███████  ██  ██   ██  ██     ██     ██      
# ██ █ ██  ██  ██   ██  ██     ██     ████    
# ██   ██  ██  ██   ██  ██     ██     ██      
# ██   ██  ██  ██     ███      ██     ██      
# ██   ██   ████      ██      ████    ██████  

def update_upw_dow_mov(G,rule, global_var):

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
    a_Ni = rule["a_Ni"]
    a_B = rule["a_B"]
    BG_thr = rule["BG_thr"]
    Ni_thr = rule["Ni_thr"]
    Bi_thr = rule["Bi_thr"]
    max_behavior = rule["max_behavior"]

    health_status = global_var.health_status
    behaviors = global_var.behavior

    # check if there are infected nodes, if not the health update will be skipped
    N_infected = np.sum(health_status == 2)
    if N_infected == 0:
        N_protecting = np.sum(np.array(behaviors == True))   #check if there are protecting nodes, if not, everything will be skipped
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

        global_var.health_status = update_health_SIR(N_infected, global_var)

        global_var.I_peak = max(global_var.I_peak  ,  float(np.mean(global_var.health_status==2)))

    N_infected = np.sum(np.array(global_var.health_status==2))
    global_var.first_tick = False
    #------------------------------------------------------------------------------------------------------------------------------#
    #third: calculate update of protection probabilities

    probability = global_var.functions.calc_protection_probability(behaviors, N_infected, global_var,
                                                                   a_B = a_B, a_Ni = a_Ni, mu = mu,
                                                                   theta = Ni_thr, Bi_thr = Bi_thr, BG_thr = BG_thr)
    
    #------------------------------------------------------------------------------------------------------------------------------#
    #fourth: update the betas and awarenesses

    global_var.behavior, global_var.beta = update_beta(probability, N_infected, max_behavior, beta0, global_var.N_nodes_H)

  
# ██  ██   ██████     ██     ██  ██    ████     ████     ████    ████     ██████  
# ██  ██   ██        ████    ██  ██     ██     ██  ██     ██     ██ ██    ██      
# ██  ██   ██       ██  ██   ██  ██     ██     ██         ██     ██  ██   ██      
# ██████   ████     ██████   ██  ██     ██      ████      ██     ██  ██   ████    
# ██  ██   ██       ██  ██   ██  ██     ██         ██     ██     ██  ██   ██      
# ██  ██   ██       ██  ██     ███      ██     ██  ██     ██     ██ ██    ██      
# ██  ██   ██████   ██  ██     ██      ████     ████     ████    ████     ██████  


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
    N_infected = np.sum(np.array(global_var.health_status))
    if N_infected == 0:
        N_protecting = np.sum(np.array(global_var.behavior==1))   #check if there are protecting nodes, if not, everything will be skipped
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
        I2R = global_var.I2R
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(len(g_h.vs)))    
        #each node's neighbors
        all_h_neighbors = global_var.h_neighbors
        #the number of infected neighbors for each node
        infected_nghbrs = np.array([np.sum(health_status[i] == 2) for i in all_h_neighbors])
        #the infection probability for each node, as if it were susceptible
        infection_probability = 1 - np.power (1 - beta,infected_nghbrs)

        health_status = ( health_status                                        +   #previous health_status
                          (rr < I2R)                    * (health_status == 2) +   #+1 if inf. agent recovers
                          (rr < infection_probability)  * (health_status == 1)     #+1 if susc. agent is infected
                        )

        global_var.health_status = health_status

        g_h.vs["health_status"] = health_status.tolist()
        #writing the health status to nodes, because I still use that when saving

        global_var.I_peak = max(global_var.I_peak  ,  float(np.mean(health_status==2)))

def update_downward_Heav(G,rule, global_var):

    """this simulates an SIR model, an assimilation model desgined to compare to an upward threshold

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
        I2R = global_var.I2R
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(len(g_h.vs)))    
        #each node's neighbors
        all_h_neighbors = global_var.h_neighbors
        #the number of infected neighbors for each node
        infected_nghbrs = np.array([np.sum(health_status[i] == 2) for i in all_h_neighbors])
        #the infection probability for each node, as if it were susceptible
        infection_probability = 1 - np.power (1 - beta,infected_nghbrs)

        health_status = ( health_status                                        +   #previous health_status
                          (rr < I2R)                    * (health_status == 2) +   #+1 if inf. agent recovers
                          (rr < infection_probability)  * (health_status == 1)     #+1 if susc. agent is infected
                        )

        global_var.health_status = health_status

        g_h.vs["health_status"] = health_status.tolist()
        #writing the health status to nodes, because I still use that when saving

        global_var.I_peak = max(global_var.I_peak  ,  float(np.mean(health_status==2)))

#  ████    ██  ██    ████    ██████  
#   ██     ███ ██     ██       ██    
#   ██     ██████     ██       ██    
#   ██     ██████     ██       ██    
#   ██     ██ ███     ██       ██    
#   ██     ██  ██     ██       ██    
#  ████    ██  ██    ████      ██    


def init_mix_3_populations(P_dyn, G, global_var):

    global_var.first_tick = True
    global_var.P_dyn = P_dyn

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    global_var.b_neighbors = [np.array(G[bl].neighbors(i)) for i in range(G[bl].vcount())]
    global_var.B_neighbor_indexing = np.array(global_var.b_neighbors) if is_regular(G[bl]) else None
    global_var.h_neighbors = [np.array(G[hl].neighbors(i)) for i in range(G[hl].vcount())]
    global_var.functions.row_mean_B = [row_mean_ragged,row_mean_fast][is_regular(G[bl])]
    global_var.functions.row_mean_H = [row_mean_ragged,row_mean_fast][is_regular(G[hl])]
    global_var.functions.calc_protection_probability = [calc_protection_probability_irregular,calc_protection_probability_regular][is_regular(G[bl])]
    global_var.N_nodes_B = G[bl].vcount()
    global_var.N_nodes_H = G[hl].vcount()

    global_var.I2R  = P_dyn["HEALTH"]["I2R"]                                                       # for each node sets the gamma

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    global_var.health_status = np.array(G[hl].vs["health_status"])
    G[bl].vs["beta"] = [P_dyn["HEALTH"]["beta0"]]*len(G[bl].vs) #  list(np.full( shape=len(G[bl].vs), fill_value = beta0))
    G[bl].vs["behavior"] = np.zeros( shape=len(G[bl].vs) )
    G[bl].vs["next_beta"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])
    global_var.I_peak = 0

    #calculating the size of each fraction in the population
    share_herders = P_dyn["BEHAVIOR"]["IC"].get("share_herders",0)
    share_contrarians = P_dyn["BEHAVIOR"]["IC"].get("share_contrarians",0)
    share_static = P_dyn["BEHAVIOR"]["IC"].get("share_static",0)

    #also filling up to 1, if possible
    if abs(share_contrarians + share_herders + share_static - 1) > 1e-6:
        fixable = False

        if P_dyn["BEHAVIOR"]["IC"].get("undefined") == "static":
            if (share_contrarians + share_herders) < 1:            
                share_static = 1 - (share_contrarians + share_herders)
                fixable = True
        elif P_dyn["BEHAVIOR"]["IC"].get("undefined") == "contrarians":
            if (share_static + share_herders) < 1:            
                share_contrarians = 1 - (share_static + share_herders)
                fixable = True
        elif P_dyn["BEHAVIOR"]["IC"].get("undefined") == "herders":
            if (share_static + share_contrarians) < 1:            
                share_herders = 1 - (share_static + share_contrarians)
                fixable = True

        if not fixable:
            print("fractions of agents are not adding up to one.")
            exit()

    #assigning all vertices an IRF group.
    remaining_idcs = np.array(range(G[bl].vcount()))
    aligner_idcs = np.random.choice(G[bl].vcount(),size=int(round(share_herders*G[bl].vcount())),replace=False)
    remaining_idcs = np.setdiff1d(remaining_idcs, aligner_idcs)
    if remaining_idcs.shape[0] != 0 and share_contrarians != 0:
        if share_static == 0:
            #share_contrarians is >0 , share_static is zero.
            #Then, make sure that no statics are assigned to the population.
            contrarian_idcs = remaining_idcs
        else:
            contrarian_idcs = np.random.choice(remaining_idcs,size=min(len(remaining_idcs),int(round(share_contrarians*G[bl].vcount()))),replace=False)
        remaining_idcs = np.setdiff1d(remaining_idcs, contrarian_idcs)
    else:
        contrarian_idcs = np.array([],dtype=bool)

    #setting up the attribute that the Metropolis algorithm will distribute
    global_var.IRF_group = np.zeros(shape=G[bl].vcount())
    global_var.IRF_group[contrarian_idcs] = 1
    global_var.IRF_group[remaining_idcs] = 2

    #if homophily flag is True and the population is not homogeneous, RE-distribute the attribute "IRF_group"
    if P_dyn["BEHAVIOR"]["IC"]["homophily"]["Flag"] == True and share_herders != 1 and share_contrarians != 1 and share_static != 1:  #not doing calculations in homogeneous populations
        set_initial_condition(P_dyn["BEHAVIOR"]["IC"],
                              attribute_vector = global_var.IRF_group,
                              g = G[bl],
                              global_var = global_var)
    
    #after potentially homophilously distriubting the attribute, write the indices to global_var
    global_var.aligner_idcs = np.where(global_var.IRF_group == 0)[0]
    global_var.contrarian_idcs = np.where(global_var.IRF_group == 1)[0]
    global_var.static_idcs = np.where(global_var.IRF_group == 2)[0]
        
    global_var.static_probability = P_dyn["BEHAVIOR"]["static_probability"]
    global_var.probability = np.where(global_var.IRF_group == 2, global_var.static_probability, 0)
    global_var.behavior = np.array(G[bl].vs["behavior"])

    rule  = {
        'func': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'mu' : P_dyn["BEHAVIOR"]["mu"],
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'a_B' : P_dyn["BEHAVIOR"]["a_B"],    
        'a_Ni' : P_dyn["BEHAVIOR"]["a_Ni"],
        'BG_thr' : P_dyn["BEHAVIOR"].get("BG_thr",P_dyn["BEHAVIOR"]["pn_thr"]),
        'Bi_thr' : P_dyn["BEHAVIOR"]["Bi_thr"],
        'Ni_thr' : P_dyn["BEHAVIOR"]["Ni_thr"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }

    if P_dyn["BEHAVIOR"]["IC"].get("equilibrium_flag",False):
        equil_steps = P_dyn["BEHAVIOR"]["IC"].get("equil_steps",10)

        for eq_step in range(equil_steps):
            update_upw_dow_health_fixed(G,rule,global_var)
            
    G[bl].vs["behavior"] = global_var.behavior.tolist()  
    global_var.beta_override =  ((1-P_dyn["BEHAVIOR"]["max_behavior"]*global_var.behavior)*P_dyn["HEALTH"]["beta0"])
    global_var.actual_health_behavior = global_var.behavior

    return rule


def init_vignette(P_dyn, G, global_var):

    global_var.first_tick = True
    global_var.P_dyn = P_dyn

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    global_var.b_neighbors = [np.array(G[bl].neighbors(i)) for i in range(G[bl].vcount())]
    global_var.B_neighbor_indexing = np.array(global_var.b_neighbors) if is_regular(G[bl]) else None
    global_var.h_neighbors = [np.array(G[hl].neighbors(i)) for i in range(G[hl].vcount())]
    global_var.functions.row_mean_B = [row_mean_ragged,row_mean_fast][is_regular(G[bl])]
    global_var.functions.row_mean_H = [row_mean_ragged,row_mean_fast][is_regular(G[hl])]
    global_var.functions.calc_protection_probability = [calc_protection_probability_irregular,calc_protection_probability_regular][is_regular(G[bl])]
    global_var.N_nodes_B = G[bl].vcount()
    global_var.N_nodes_H = G[hl].vcount()

    global_var.I2R  = P_dyn["HEALTH"]["I2R"]                                                       # for each node sets the gamma

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    global_var.health_status = np.array(G[hl].vs["health_status"])
    G[bl].vs["beta"] = [P_dyn["HEALTH"]["beta0"]]*len(G[bl].vs) #  list(np.full( shape=len(G[bl].vs), fill_value = beta0))
    G[bl].vs["behavior"] = np.zeros( shape=len(G[bl].vs) )
    G[bl].vs["next_beta"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])
    global_var.I_peak = 0

    #calculating the size of each fraction in the population
    size_herders =     int(round(P_dyn["BEHAVIOR"]["IC"].get("share_herders",0)  *  G[bl].vcount()))
    size_contrarians = int(round(P_dyn["BEHAVIOR"]["IC"].get("share_contrarians",0)  *  G[bl].vcount()))
    size_always =      int(round(P_dyn["BEHAVIOR"]["IC"].get("share_always",0)  *  G[bl].vcount()))
    size_never =       int(round(P_dyn["BEHAVIOR"]["IC"].get("share_never",0)  *  G[bl].vcount()))
    size_medium =      int(round(P_dyn["BEHAVIOR"]["IC"].get("share_static",0)  *  G[bl].vcount()))

    size_accounted_for = size_contrarians + size_herders + size_medium + size_always + size_never
    
    #filling up to 100%, if possible
    if size_accounted_for < G[bl].vcount():
        size_unaccounted_for = G[bl].vcount() - size_accounted_for

        if P_dyn["BEHAVIOR"]["IC"].get("undefined") == "static":
            size_medium += size_unaccounted_for
        elif P_dyn["BEHAVIOR"]["IC"].get("undefined") == "contrarians":      
            size_contrarians += size_unaccounted_for
        elif P_dyn["BEHAVIOR"]["IC"].get("undefined") == "herders":         
            size_herders += size_unaccounted_for
        elif P_dyn["BEHAVIOR"]["IC"].get("undefined") == "never":       
            size_never += size_unaccounted_for
        elif P_dyn["BEHAVIOR"]["IC"].get("undefined") == "always":       
            size_always += size_unaccounted_for

    elif size_accounted_for > G[bl].vcount():
        print("fractions of agents are not adding up to one.")
        exit()

    #assigning all vertices an IRF group.
    remaining_idcs = np.array(range(G[bl].vcount()))
    aligner_idcs = np.random.choice(G[bl].vcount(),size= size_herders,replace=False)
    remaining_idcs = np.setdiff1d(remaining_idcs, aligner_idcs)

    contrarian_idcs = np.random.choice(remaining_idcs,size=size_contrarians,replace=False)
    remaining_idcs = np.setdiff1d(remaining_idcs, contrarian_idcs)

    medium_protection_idcs = np.random.choice(remaining_idcs,size=size_medium,replace=False)
    remaining_idcs = np.setdiff1d(remaining_idcs, medium_protection_idcs)

    always_protection_idcs = np.random.choice(remaining_idcs,size=size_always,replace=False)
    remaining_idcs = np.setdiff1d(remaining_idcs, always_protection_idcs)

    never_protection_idcs = remaining_idcs

    #setting up the attribute that the Metropolis algorithm will distribute
    global_var.IRF_group = np.zeros(shape=G[bl].vcount())
    global_var.IRF_group[contrarian_idcs] = 1
    global_var.IRF_group[medium_protection_idcs] = 2
    global_var.IRF_group[always_protection_idcs] = 3
    global_var.IRF_group[never_protection_idcs] = 4

    #if homophily flag is True and the population is not homogeneous, RE-distribute the attribute "IRF_group"
    if P_dyn["BEHAVIOR"]["IC"]["homophily"]["Flag"] == True:
        if size_herders != G[bl].vcount() and size_contrarians != G[bl].vcount() and size_medium != G[bl].vcount() and size_always != G[bl].vcount() and size_never != G[bl].vcount():
            #not doing calculations in homogeneous populations

            set_initial_condition(P_dyn["BEHAVIOR"]["IC"],
                                  attribute_vector = global_var.IRF_group,
                                  g = G[bl],
                                  global_var = global_var)

    #after potentially homophilously distriubting the attribute, write the indices to global_var
    global_var.aligner_idcs = np.where(global_var.IRF_group == 0)[0]
    global_var.contrarian_idcs = np.where(global_var.IRF_group == 1)[0]
    global_var.medium_idcs = np.where(global_var.IRF_group == 2)[0]
    global_var.never_idcs = np.where(global_var.IRF_group == 4)[0]
    global_var.always_idcs = np.where(global_var.IRF_group == 3)[0]
    global_var.static_idcs = np.where(global_var.IRF_group > 1)[0]

    global_var.static_probability = np.zeros(shape = G[bl].vcount())
    global_var.static_probability[medium_protection_idcs] = np.clip(np.random.uniform(low = P_dyn["BEHAVIOR"]["static_probability_min"],
                                                                                      high = P_dyn["BEHAVIOR"]["static_probability_max"],
                                                                                      size = size_medium),
                                                                    a_min = 0, a_max = 1)
    global_var.static_probability[always_protection_idcs] = 1
    #static_probability is now
    #    zero for aligners and contrarians
    #    zero for never protecters
    #    one for always protectors
    #    between 0 and 1 for medium protectors
    
    global_var.probability = global_var.static_probability
    global_var.behavior = np.array(G[bl].vs["behavior"])

    global_var.static_probability = global_var.static_probability[global_var.static_idcs]
    #filtering out all aligners and contrarians

    rule  = {
        'func': P_dyn["func"],
        'hl': hl,
        'bl': bl,
        'mu' : P_dyn["BEHAVIOR"]["mu"],
        'beta0' : P_dyn["HEALTH"]["beta0"],
        'a_B' : P_dyn["BEHAVIOR"]["a_B"],    
        'a_Ni' : P_dyn["BEHAVIOR"]["a_Ni"],
        'BG_thr' : P_dyn["BEHAVIOR"].get("BG_thr",P_dyn["BEHAVIOR"]["pn_thr"]),
        'Bi_thr' : P_dyn["BEHAVIOR"]["Bi_thr"],
        'Ni_thr' : P_dyn["BEHAVIOR"]["Ni_thr"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }

    if P_dyn["BEHAVIOR"]["IC"].get("equilibrium_flag",False):
        equil_steps = P_dyn["BEHAVIOR"]["IC"].get("equil_steps",10)

        for eq_step in range(equil_steps):
            update_upw_dow_health_fixed(G,rule,global_var)
            
    G[bl].vs["behavior"] = global_var.behavior.tolist()
    global_var.actual_health_behavior = global_var.behavior

    return rule


def init_upw_dow(P_dyn, G,global_var):

    global_var.first_tick = True

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    global_var.functions.row_mean_B = [row_mean_ragged,row_mean_fast][is_regular(G[bl])]
    global_var.functions.row_mean_H = [row_mean_ragged,row_mean_fast][is_regular(G[hl])]
    global_var.functions.calc_protection_probability = [calc_protection_probability_irregular,calc_protection_probability_regular][is_regular(G[bl])]
                                                        
    for i,vertex in enumerate(G[hl].vs):
        vertex["b_neighbors"] = G[bl].neighbors(i)
        vertex["h_neighbors"] = G[hl].neighbors(i)
    global_var.b_neighbors = [np.array(G[bl].neighbors(i)) for i in range(G[bl].vcount())]
    global_var.B_neighbor_indexing = np.array(global_var.b_neighbors) if is_regular(G[bl]) else None
    global_var.h_neighbors = [np.array(G[hl].neighbors(i)) for i in range(G[hl].vcount())]
    global_var.N_nodes_B = G[bl].vcount()
    global_var.N_nodes_H = G[hl].vcount()

    global_var.I2R = P_dyn["HEALTH"]["I2R"]
    #just P_dyn["HEALTH"]["I2R"]; this line is complicated only for backwards compatibility

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    global_var.health_status = np.array(G[hl].vs["health_status"])
    global_var.beta = [P_dyn["HEALTH"]["beta0"]]*len(G[bl].vs) #  list(np.full( shape=len(G[bl].vs), fill_value = beta0))
    global_var.behavior = np.full( shape=len(G[bl].vs), fill_value = 0).tolist()
    G[bl].vs["next_beta"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])
    global_var.I_peak = 0

    if "UPW" in P_dyn["func"]:
        global_var.aligner_idcs = range(global_var.N_nodes_B)
        global_var.contrarian_idcs = []
        global_var.static_idcs = []
        global_var.IRF_group = np.zeros(shape = G[bl].vcount())

        G[bl].vs["herder"] = True
    else:
        G[bl].vs["herder"] = False
        global_var.aligner_idcs = []
        global_var.contrarian_idcs = range(global_var.N_nodes_B)
        global_var.static_idcs = []
        global_var.IRF_group = np.ones(shape = G[bl].vcount())
        
    G[bl].vs[global_var.aligner_idcs]["IRF_group"] = 0
    G[bl].vs[global_var.contrarian_idcs]["IRF_group"] = 1
    G[bl].vs[global_var.static_idcs]["IRF_group"] = 2

    global_var.static_probability = P_dyn["BEHAVIOR"].get("static_probability",None)

    if P_dyn["BEHAVIOR"]["IC"].get("equilibrium_flag",False):
        g_b = G[bl]
        a_B = P_dyn["BEHAVIOR"]["a_B"]
        mu = P_dyn["BEHAVIOR"]["mu"]
        BG_thr = P_dyn["BEHAVIOR"].get("BG_thr",P_dyn["BEHAVIOR"]["pn_thr"])
        N_infected = P_dyn["HEALTH"]["IC"]["N_pat_zero"]
        #g_h = G[hl]
        Bi_thr = P_dyn["BEHAVIOR"]["Bi_thr"]
        a_Ni = P_dyn["BEHAVIOR"]["a_Ni"]
        Ni_thr = P_dyn["BEHAVIOR"]["Ni_thr"]
        max_behavior = P_dyn["BEHAVIOR"]["max_behavior"]
        beta0 = P_dyn["HEALTH"]["beta0"]
        equil_steps = P_dyn["BEHAVIOR"]["IC"].get("equil_steps",10)

        behavior = np.array(global_var.behavior)

        for eq_step in range(equil_steps):

            probability = global_var.functions.calc_protection_probability(behavior, N_infected, global_var,
                                                      a_B = a_B, a_Ni = a_Ni, mu = mu,
                                                      theta = Ni_thr, Bi_thr = Bi_thr, BG_thr = BG_thr)
            behavior, beta = update_beta(probability, N_infected, max_behavior, beta0, g_b.vcount())

        global_var.behavior = behavior
        global_var.beta = beta
    global_var.beta_override = beta

    if "UPW" in P_dyn["func"]:
        IRF_direction = "UPW"
    elif "DOW" in P_dyn["func"]:
        IRF_direction = "DOW"
    else:
        IRF_direction = None

    rule  = {
        'func': P_dyn["func"],
        'IRF_direction' : IRF_direction,
        'hl': hl,
        'bl': bl,
        'mu' : P_dyn["BEHAVIOR"]["mu"],
        'beta0' : P_dyn["HEALTH"]["beta0"],
        "a_B" : P_dyn["BEHAVIOR"]["a_B"],
        'a_Ni' : P_dyn["BEHAVIOR"]["a_Ni"],
        'BG_thr' : P_dyn["BEHAVIOR"].get("BG_thr",P_dyn["BEHAVIOR"]["pn_thr"]),
        'Bi_thr' : P_dyn["BEHAVIOR"]["Bi_thr"],
        'Ni_thr' : P_dyn["BEHAVIOR"]["Ni_thr"],
        'max_behavior' :  P_dyn["BEHAVIOR"]["max_behavior"]
        }
    return rule

def init_Heav(P_dyn, G,global_var):

    global_var.first_tick = True

    hl = P_dyn["HEALTH"]["layer"]       # layer where the the health status is imprinted
    bl = P_dyn["BEHAVIOR"]["layer"]     # layer where the behavior is imprinted

    global_var.I2R = P_dyn["HEALTH"]["I2R"]

    # for each node sets the initial condition
    set_disease_initial_condition(P_dyn["HEALTH"]["IC"], "health_status", G[hl])
    G[hl].vs["next_health"] = G[hl].vs["health_status"]
    G[bl].vs["beta"] = [P_dyn["HEALTH"]["beta0"]]*len(G[bl].vs) #  list(np.full( shape=len(G[bl].vs), fill_value = beta0))
    G[bl].vs["behavior"] = np.full( shape=len(G[bl].vs), fill_value = 0)
    G[bl].vs["next_beta"] = np.full( shape=len(G[bl].vs), fill_value = P_dyn["HEALTH"]["beta0"])
    global_var.I_peak = 0

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

def init_model(update_fct_dict, init_fct_dict):
    #main upwards and downwards IRFs
    init_fct_dict["UPW"] = init_upw_dow
    init_fct_dict["DOW"] = init_upw_dow
    update_fct_dict["UPW"] = update_upw_dow
    update_fct_dict["DOW"] = update_upw_dow

    #upwards and downwards SIS functions
    init_fct_dict["UPW_SIS"] = init_upw_dow
    init_fct_dict["DOW_SIS"] = init_upw_dow
    update_fct_dict["UPW_SIS"] = update_upw_dow_SIS
    update_fct_dict["DOW_SIS"] = update_upw_dow_SIS

    #functions for movies
    init_fct_dict["UPW_MOV"] = init_upw_dow
    init_fct_dict["DOW_MOV"] = init_upw_dow
    init_fct_dict["mix_3_MOV"] = init_mix_3_populations
    update_fct_dict["UPW_MOV"] = update_upw_dow_mov
    update_fct_dict["DOW_MOV"] = update_upw_dow_mov    
    update_fct_dict["mix_3_MOV"] = update_upw_dow_mov

    #step function IRFs
    init_fct_dict["UPW_Heav"] = init_Heav
    init_fct_dict["DOW_Heav"] = init_Heav
    update_fct_dict["UPW_Heav"] = update_upward_Heav  
    update_fct_dict["DOW_Heav"] = update_downward_Heav
    
    #mixed populations
    init_fct_dict["mix_3"] = init_mix_3_populations
    update_fct_dict["mix_3"] = update_upw_dow
    init_fct_dict["vignette"] = init_vignette
    update_fct_dict["vignette"] = update_upw_dow






