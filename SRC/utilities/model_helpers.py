import numpy as np

# ██  ██   ██████   ██       █████    ██████   █████             ██████   ██  ██   ██  ██    ████    ██████    ████     ████    ██  ██    ████   
# ██  ██   ██       ██       ██  ██   ██       ██  ██            ██       ██  ██   ███ ██   ██  ██     ██       ██     ██  ██   ███ ██   ██  ██  
# ██  ██   ██       ██       ██  ██   ██       ██  ██            ██       ██  ██   ██████   ██         ██       ██    ██    ██  ██████   ██      
# ██████   ████     ██       █████    ████     █████             ████     ██  ██   ██████   ██         ██       ██    ██    ██  ██████    ████   
# ██  ██   ██       ██       ██       ██       ████              ██       ██  ██   ██ ███   ██         ██       ██    ██    ██  ██ ███       ██  
# ██  ██   ██       ██       ██       ██       ██ ██             ██       ██  ██   ██  ██   ██  ██     ██       ██     ██  ██   ██  ██   ██  ██  
# ██  ██   ██████   ██████   ██       ██████   ██  ██            ██        ████    ██  ██    ████      ██      ████     ████    ██  ██    ████  


def calc_protection_probability_regular(behaviors, N_infected, global_var,
                                a_B = 6, a_Ni = 25, mu = 0.75,
                                theta = 0.04, Bi_thr = 0.5, BG_thr = 0.5):
    
    protecting_nghbrs = np.mean(behaviors[global_var.B_neighbor_indexing], axis = 1)

    exponent = np.zeros(global_var.N_nodes_B)

    exponent_upw = (- a_B * mu     * (protecting_nghbrs       - BG_thr)
                    - a_B * (1-mu) * (behaviors      - Bi_thr) 
                    - a_Ni*          (N_infected/global_var.N_nodes_H - theta))
    
    exponent_dow = (+ a_B * mu     * (protecting_nghbrs       - BG_thr)
                    - a_B * (1-mu) * (behaviors      - Bi_thr)
                    - a_Ni*          (N_infected/global_var.N_nodes_H - theta))
    
    exponent[global_var.aligner_idcs] = exponent_upw[global_var.aligner_idcs]
    exponent[global_var.contrarian_idcs] = exponent_dow[global_var.contrarian_idcs]
    #implicitly having the exponent for the remaining indices =0
    #I overwrite the result of this exponent two lines beneath in the probability vector

    probability = 1 / (1 + np.exp( exponent ) )
    probability[global_var.static_idcs] = global_var.static_probability

    return probability


def calc_protection_probability_irregular(behaviors, N_infected, global_var,
                                a_B = 6, a_Ni = 25, mu = 0.75,
                                theta = 0.04, Bi_thr = 0.5, BG_thr = 0.5):

    protecting_nghbrs = np.array([np.mean(behaviors[i]) for i in global_var.b_neighbors])

    exponent = np.zeros(global_var.N_nodes_B)

    exponent_upw = (- a_B * mu     * (protecting_nghbrs       - BG_thr)
                    - a_B * (1-mu) * (behaviors      - Bi_thr) 
                    - a_Ni*          (N_infected/global_var.N_nodes_H - theta))
    
    exponent_dow = (+ a_B * mu     * (protecting_nghbrs       - BG_thr)
                    - a_B * (1-mu) * (behaviors      - Bi_thr)
                    - a_Ni*          (N_infected/global_var.N_nodes_H - theta))
    
    exponent[global_var.aligner_idcs] = exponent_upw[global_var.aligner_idcs]
    exponent[global_var.contrarian_idcs] = exponent_dow[global_var.contrarian_idcs]
    #implicitly having the exponent for the remaining indices =0
    #I overwrite the result of this exponent two lines beneath in the probability vector

    probability = 1 / (1 + np.exp( exponent ) )
    probability[global_var.static_idcs] = global_var.static_probability

    return probability

def update_behavior(probability, N_infected, N_nodes):
    rr1 = np.random.uniform(low=0, high=1, size=N_nodes)
    
    behavior = ( (probability > rr1) * (N_infected > 0) )
    #                       dice throw          if there are infected  agents
    #agents will stop protecting if no one is infected
    #this is an assumption/approximation that doesn't change any interesting predictions but speeds up simulations dramatically

    return behavior

def update_beta(probability, N_infected, max_behavior, beta0, N_nodes):
    rr1 = np.random.uniform(low=0, high=1, size=N_nodes)
    
    behavior = ( (probability > rr1) * (N_infected > 0) )
    #                       dice throw          if there are infected  agents
    #agents will stop protecting if no one is infected
    #this is an assumption/approximation that doesn't change any interesting predictions but speeds up simulations dramatically
    beta = (1-max_behavior*behavior) * beta0

    return behavior, beta

def update_health_SIR(N_infected, global_var):
    health_status = global_var.health_status

    if(N_infected>0):
        I2R = global_var.I2R
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(global_var.N_nodes_H))

        #the number of infected neighbors for each node
        infected_nghbrs = np.array([np.sum(health_status[i] == 2) for i in global_var.h_neighbors])

        #the infection probability for each node, as if it were susceptible
        infection_probability = 1 - np.power (1 - global_var.beta,infected_nghbrs)

        return (   health_status                                          #previous health_status
                 + (rr < I2R)                    * (health_status == 2)   #+1 if inf. agent recovers
                 + (rr < infection_probability)  * (health_status == 1)   #+1 if susc. agent is infected
               )
    
    return global_var.health_status

def update_health_SIS(N_infected, global_var):
    health_status = global_var.health_status

    if(N_infected>0):
        I2R = global_var.I2R
        
        # generate a random number for each node
        rr=np.random.uniform(low=0, high=1, size=(global_var.N_nodes_H))

        #the number of infected neighbors for each node
        infected_nghbrs = np.array([np.sum(health_status[i] == 2) for i in global_var.h_neighbors])

        #the infection probability for each node, as if it were susceptible
        infection_probability = 1 - np.power (1 - global_var.beta,infected_nghbrs)

        return (   health_status                                          #previous health_status
                 - (rr < I2R)                    * (health_status == 2)   #-1 if infected agent recovers and is susceptible again
                 + (rr < infection_probability)  * (health_status == 1)   #+1 if susc. agent is infected
               )
    
    return global_var.health_status