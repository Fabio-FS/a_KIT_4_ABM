#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import igraph as ig
import random as rn
import pathlib
#from numba import jit
import warnings

from Parameter_Importer import *
from Utility_Functions  import *
from Clustering import *
    
    
    
    
    
    

    
#########################################
###The Big 6#############################
#########################################

#global infected

def Global_fraction_disease(g,n,h_s=2):   #standard target is 2 for Infected
    avg_D = np.sum([1 for v in g.vs if v["health_status"] == h_s]) / n
    return avg_D


def Global_f_behavior(g,PS,f=np.average):   #standard function is np.average
    #avg_D = np.sum([1 for v in g.vs if v["health_status"] == h_s]) / PS.n
    GFB = f([1 for v in g.vs if v["health_status"] == h_s])
    return GFB


def local_fraction_disease(g,PS,h_s=2, default_value=0):
    #with warnings.catch_warnings():
    #    warnings.simplefilter("ignore", category=RuntimeWarning)

    #takes the average behavior for neighbors IF TOO SLOW I HAVE SAVED THE LIST OF NEIGHBORS FOR EACH NODE IN PD.neighbors
    g.vs['L_avg_D'] = np.zeros(PS.n)

    for v in g.vs:
        if(v['L_dis_neighbors']>0):
            v['L_avg_D'] = len(my_filtering_target(g, v, "dis_neighbors", "health_status", h_s))/v['L_dis_neighbors']
        else:
            v["L_avg_D"] = default_value

def local_f_behavior(g,PS,PB, f =np.average, default_value=0.5):    
            
        #takes the average behavior for neighbors IF TOO SLOW I HAVE SAVED THE LIST OF NEIGHBORS FOR EACH NODE IN PD.neighbors
        g.vs['L_f_B'] = np.zeros(PS.n)
        
        for v in g.vs:
            if(v["L_beh_neighbors"]>0):
                v['L_f_B'] = f([v2['behavior_status'] for v2 in g.vs[v['beh_neighbors']]])
            else:
                v["L_f_B"] = default_value
                
def local_f_disease(g,PS,PB, f =np.average, default_value=0.5):    
            
        #takes the average behavior for neighbors IF TOO SLOW I HAVE SAVED THE LIST OF NEIGHBORS FOR EACH NODE IN PD.neighbors
        g.vs['L_f_D'] = np.zeros(PS.n)
        
        for v in g.vs:
            if(v["L_beh_neighbors"]>0):
                v['L_f_D'] = f([v2['health_status'] for v2 in g.vs[v['dis_neighbors']]])
            else:
                v["L_f_D"] = default_value
    

def single_target_status(g,PS,PX,attribute_name,full_name,default_value=0):
    if((attribute_name == 'beh')|(attribute_name == 'dis')):
        name='single_target_' + attribute_name
        #print(name)
        name_type = attribute_name + '_neighbors'
        for v in g.vs:
            if(v['L_' + name_type]>0):
                v[name] = g.vs[rn.choice(v[name_type])][full_name]
            else:
                v[name] = default_value
    else:
        print("attribute_name=" + attribute_name)
        input('attribute name wrong in single target status! only beh or dis expected!')
    
#########################################################
#                    BASIC     IRF                      #
#########################################################


def evolve_Average(g,PS,PD,PB):
    
    # in this scenario agent i copies the behavior of other agents if there are enough of them behaving in that way + random    
    g.vs['new_state_B']=np.zeros(PS.n)
    g.vs['L_f_B']=np.zeros(PS.n)

    single_target_status(g,PS,PB,'beh','behavior_status',default_value = 0.5)   #updates v['single_target_beh']
    #local_f_behavior(g,PS,PB) #  it updates v['L_f_B']

    epsilon = 0.1
   
    ## the update rule is: B_i(t+1) = K * (B_i(t)) + (1-K) * (B_j-B_i)
    for v in g.vs:
        v['new_state_B'] = v['behavior_status'] + 0.1*(v['single_target_beh']-v['behavior_status'])


#bounded confidence model
def evolve_BCM(g,PS,PD,PB):
    # For the Bounded confidence model see Deffuant et al. 2000, Hegselmann & Krause 2002
    g.vs['new_state_B']=np.zeros(PS.n)
    g.vs['L_f_B']=np.zeros(PS.n)

    single_target_status(g,PS,PB,'beh','behavior_status',default_value = 0.5)   #updates v['single_target_beh']

    epsilon = 0.1
    mu = 0.1
    for v in g.vs:
        i = v['behavior_status']
        j = v['single_target_beh']
        # epsilon is the confidence level (Hegselmann & Krause 2002)
        # mu is the rate of opinion convergence 0<mu<0.5
        if(np.abs(i-j)>epsilon):
            v['new_state_B'] = v['behavior_status'] + mu*(j-i)
        else:
            v['new_state_B'] = v['behavior_status']
            

def evolve_Negative_Influence(g,PS,PD,PB):
    #Jager & Amblard 2005

    g.vs['new_state_B']=np.zeros(PS.n)
    g.vs['L_f_B']=np.zeros(PS.n)

    single_target_status(g,PS,PB,'beh','behavior_status',default_value = 0.5)   #updates v['single_target_beh']

    mu = 0.1
    u = 0.3
    for v in g.vs:
        i = v['behavior_status']
        j = v['single_target_beh']
        # mu is the rate of opinion convergence 0<mu<0.5
        if(np.abs(i-j)>u):
            v['new_state_B'] = v['behavior_status'] + mu*(j-i)
        else:
            v['new_state_B'] = v['behavior_status'] - mu*(j-i)
        v['new_state_B'] = min(max(v['new_state_B'],0),1)
            





def evolve_Reinforcement(g,PS,PD,PB):
    # for simplicity I keep it between -1 and 1, but it can be easily rescaled between 0 and 1.
 
    g.vs['new_state_B']=np.zeros(PS.n)
    
    single_target_status(g,PS,PB,'beh','behavior_status')   # updates v['single_target_beh']
                                       
    delta = 0.1
    
    for v in g.vs:
        
        i = 2*v['behavior_status']-1
        j = 2*v['single_target_beh']-1

        # if they have the same predisposition, i is reinforced in its own direction. (respect to the null opinion 0.5)

        v['new_state_B'] = i +delta*(np.heaviside(np.sign(i)*np.sign(j),0))*np.sign(i)
        v['new_state_B'] = min(max(v['new_state_B'],0),1)

# Forgot about this
def evolve_Repulsive(g,PS,PD,PB):
    # opinion between 0 and 1. If the difference is bigger than threshold, then the agent's opinion moves away from target, otherwise it moves towards it.

    Threshold = 0.25
    delta = 0.05
    
    g.vs['new_state_B']=np.zeros(PS.n)
    
    single_target_status(g,PS,PB,'beh','behavior_status')   # updates v['single_target_beh']


    for v in g.vs:
        i = v['behavior_status']
        j = v['single_target_beh']
        variation = delta*np.sign(j-i)*np.sign(Threshold-abs(i-j))
        
        v['new_state_B'] = i + variation
        v['new_state_B'] = min(max(v['new_state_B'],0),1)


#########################################################
#   INFLUENCE RESPONSE FUNCTIONS FROM LITTERATURE       #
#########################################################

    
    
    
def evolve_Xu2022(g,PS,PD,PB):
    """
    In this scenario each agent has a behavior (level of awareness) that depends on:
    1) his connectivity k,
    2) the fraction of first neighboros infected
    3) the global fraction of infected.
    
    B = 0.5^(k*(0.2*L+0.8*I))
    
    """
    
    g.vs['new_state_B']=np.zeros(PS.n)
    avg_D = Global_fraction_disease(g,PS)
    local_fraction_disease(g,PS) #update v['L_avg_D']
    
    
    p = 0.2 # if interested put it in the PB, and access it from there.
    
    for v in g.vs:
        v['new_state_B'] = np.power(0.5, v['L_dis_neighbors'] * (p*v['L_avg_D']+(1-p)*avg_D))


    #super slow and most general scenario, useful to copy what is needed and adapt.
    
    
    
def evolve_Xu2022_light(g,PS,PD,PB):
    """
    In this scenario each agent has a behavior (level of awareness) that depends on:
    1) his connectivity k,
    2) the global fraction of infected.
    
    B = 0.5^(kI)
    
    """
    
    g.vs['new_state_B']=np.zeros(PS.n)
    avg_D = Global_fraction_disease(g,PS)

    
    for v in g.vs:
        v['new_state_B'] = np.power(0.5, v['L_dis_neighbors']*avg_D)
        v['new_state_B'] = 1

    
    #super slow and most general scenario, useful to copy what is needed and adapt.
    
    
def evolve_XIA_2013(g,PS,PD,PB):
    for v in g.vs:
        v['health_status']=v["behavior_status"]
        
    ## p = 1 always game theory, p = 0, always peer pressure
    
    
    
    #g.vs['new_state_B']=np.zeros(PS.n)
    
    local_f_disease(g, PS, PB, f =np.sum, default_value=1)   # saves it in g.vs['L_f_D']
    # Here the sum over the nodes of the health status is the numbber of vaccinated
    
    rr = np.random.uniform(0,1,size = PS.n) # choice between peer pressure and gametheory
    rr2 = np.random.uniform(0,1,size = PS.n) # if game theory I need a random number

    STRATEGY = np.zeros(PS.n)
    STRATEGY[rr<PB.p]=1
    
    
    
    
    for i,v in enumerate(g.vs):
        if(STRATEGY[i] == 0):
            # use peer pressure
            if(v['L_dis_neighbors'] == 0):
                v['new_state_B'] = 1
            else:
                sq1 = np.sqrt(v["L_f_D"])
                sq2 = np.sqrt(v['L_dis_neighbors']-v["L_f_D"])
                
                DELTA = (sq1-sq2)/(sq1 + sq2)
                
                P = 1/(1+np.exp(-DELTA*PB.alpha))
                
                if(rr2[i]<P):
                    v['new_state_B'] = 1
                    v['health_status'] = 1
                else:
                    v['new_state_B'] = 0
                    v['health_status'] = 0           
        else:
            # use game theory
            #print("game theory")
            if(v['L_dis_neighbors'] == 0):
                v['new_state_B'] = 1
                v['health_status'] = 1
            else:
                lambda_hat = PD.S2I*(1-v["L_f_D"]/v["L_dis_neighbors"])
                if(PB.r<lambda_hat):
                    v['new_state_B'] = 1
                    v['health_status'] = 1
                else:
                    v['new_state_B'] = 0
                    v['health_status'] = 0
    for v in g.vs:
        if(v["L_dis_neighbors"]>0):
            v["L_avg_D"] = v["L_f_D"]/v["L_dis_neighbors"]
        else:
            v["L_avg_D"] = 1
    #count = 0
    #for v in g.vs:
    #    if (v["new_state_B"] == 4):
    #        count = count + 1
    #input("I am counting: " + str( count) + " vaccinated")
    
def evolve_Xia_heaviside(g,PS,PD,PB):
    
    # in this scenario agent i copies the behavior of other agents if there are enough of them behaving in that way + random    
    g.vs['new_state_B']=np.zeros(PS.n)
    
    # metacode: v(i) = p_i(np.heaviside(0.6-V))+(1-p_i)*(1/(1+np.exp(-10*(V+0.5))))
    
    avg_V = Global_fraction_disease(g,PS,h_s = PB.influencing_health_status)

    for v in g.vs:
        v['new_state_B'] = v['peer_pressure']* np.heaviside(0.6-avg_V,0.5) +(1-v['peer_pressure'])*(1/(1+np.exp(-10*(avg_V-0.5))))


def evolve_PIRES(g,PS,PD,PB):
    
    # in this scenario agent i copies the behavior of other agents if there are enough of them behaving in that way + random    
    new_state = np.zeros(PS.n)
    
    local_f_behavior(g,PS,PB)               #updates the v['L_f_B']
    avg_D = Global_fraction_disease(g,PS)
    
    epsilon = 0.1 #  in the future this could be included in the PB parameter set. But for now it stays here
    for v in g.vs:
        v['new_state_B'] = v['behavior_status'] + epsilon*v['L_f_B'] + v['global_fear']*avg_D
        v['new_state_B'] = min(max(v['new_state_B'],0),1)
    



def old_evolve_Reinforcement(g,PS,PD,PB):
    # for simplicity I keep it between -1 and 1, but it can be easily rescaled between 0 and 1.
 
    g.vs['new_state_B']=np.zeros(PS.n)
    
    single_target_status(g,PS,PB,'beh','behavior_status')   # updates v['single_target_beh']
                                       
    news = 0.02
    
    # B(t+1) = p B(t) +(1-p)[H(Bi)H(Bj)-H(-Bi)H(-Bj)]  -- rescaled into x\in[-1, +1]
    for v in g.vs:
        #print("me status= " + str(2*v['behavior_status']-1) + " him status = " + str(2*v['single_target_beh']-1))
        
        i = v['behavior_status']
        j = v['single_target_beh']
        
        
        #print("I=" + str(i) + ", J= " + str(j) + ", Delta=" + str(Delta))
        
        if((np.sign(i-0.5)==-1)&(np.sign(j-0.5)==-1)):
            Delta = -0.5
        elif((np.sign(i-0.5)==+1)&(np.sign(j-0.5)==+1)):
            Delta = +0.5
        else:
            Delta = 0
        v['new_state_B'] = i + news*Delta
        v['new_state_B'] = min(max(v['new_state_B'],0),1)
                                       
 
# %%



#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import igraph as ig
import random as rn
import pathlib
#from numba import jit

from Parameter_Importer import *
from Utility_Functions  import *




## SIR where the probability of being infected depends ONLY on the behavior of the susceptible individuals


def evolve_SIRB(g,n,PD,PB):
    
    g.vs['new_state_D']=np.zeros(n)
    
    PB = np.array(g.vs["behavior_status"])    # Personal behavior. 

    g.vs['personal_beta'] = 1 - PD.S2I_l + (1-PB)*(PD.S2I_l-PD.S2I_u)
    
    #stat = np.array(g.vs['health_status'])
    
    #n_inf=[len(np.array(my_filtering_target(g,v,"dis_neighbors","health_status",2))) for v in g.vs]
    
    #stat=np.array(stat)
    N_infected = np.sum(np.array(g.vs['health_status'])==2)
    
    if(N_infected>0):
        rr=np.random.uniform(low=0, high=1, size=(n))
        for i,vertex in enumerate(g.vs):
            if vertex['health_status']==1:
                h_n=np.array(my_filtering_target(g,vertex,"dis_neighbors","health_status",2))
                
                P_infected = np.power(vertex['personal_beta'],len(h_n))
                
                if(rr[i]<P_infected):
                    vertex['new_state_D']=1
                else:
                    vertex['new_state_D']=2
                    
            elif(vertex['health_status']==2):
                if(rr[i]<PD.I2R):
                    #input("healed")
                    vertex['new_state_D']=3
                else:
                    #input("still sick")
                    vertex['new_state_D']=2
            elif vertex['health_status']==3:
                vertex['new_state_D']=3
                
            else:
                print('ERROR! This node had an health state that was unexpected!!! in evolve_SIRB')
                wait = input("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    else:
        for i,vertex in enumerate(g.vs):
            vertex['new_state_D'] = vertex['health_status']



def evolve_SIR_B_recovery(g,n,PD,PB):
    
    g.vs['new_state_D']=np.zeros(n)
    
    PBeh = np.array(g.vs["behavior_status"])    # Personal behavior. 

    g.vs['personal_recovery'] = PD.I2R_u -PBeh*(PD.I2R_u-PD.I2R_l)
    
    #stat = np.array(g.vs['health_status'])
    
    #n_inf=[len(np.array(my_filtering_target(g,v,"dis_neighbors","health_status",2))) for v in g.vs]
    
    #stat=np.array(stat)
    N_infected = np.sum(np.array(g.vs['health_status'])==2)
    
    if(N_infected>0):
        rr=np.random.uniform(low=0, high=1, size=(n))
        for i,vertex in enumerate(g.vs):
            if vertex['health_status']==1:
                h_n=np.array(my_filtering_target(g,vertex,"dis_neighbors","health_status",2))
                
                P_infected = np.power(1-PD.S2I,len(h_n))
                
                if(rr[i]<P_infected):
                    vertex['new_state_D']=1
                else:
                    vertex['new_state_D']=2
                    
            elif(vertex['health_status']==2):
                if(rr[i]<vertex['personal_recovery']):
                    #input("healed")
                    vertex['new_state_D']=3
                else:
                    #input("still sick")
                    vertex['new_state_D']=2
            elif vertex['health_status']==3:
                vertex['new_state_D']=3
                
            else:
                print('ERROR! This node had an health state that was unexpected!!! in evolve_SIR-B-RECOVERY')
                wait = input("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    else:
        for i,vertex in enumerate(g.vs):
            vertex['new_state_D'] = vertex['health_status']


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import igraph as ig
import random as rn
import pathlib
import json

from Behavior_Functions import Global_fraction_disease
from Utility_Functions import *
#from Save_Functions import *
#from numba import jit


#&(str(getattr(self, i))[0] != "<")
## Initialize:
class Parameters:
    def _print(self):
        print("print attributes and values")
        for i in dir(self):
            if((i[0] != "_")):
                print(str(i) + " = " + str(getattr(self, i)))
    def _print_names(self):
        for i in dir(self):
            if((i[0] != "_")):
                print(str(i))
    def _list_names(self):
        res = []
        for i in dir(self):
            if((i[0] != "_")):
                res.append(str(i))
        return res
class SAVE():
    def __init__(self,dictionary):
        # import all the parameters provided, and create the appropriate attributes
        for key, value in dictionary.items():
            setattr(self, key, value)
            
        #self._print()
        #input (" " )
        if(self.TYPE == "BEH"):
            # name of the attribute,
            if(hasattr(self, "name_attribute") == False):
                self.name_attribute = "behavior_status"

            # value of the attribute, (for example in SIR, if you only want I set value_attribute = 1)    
            if(hasattr(self, "value_attribute") == False):
                self.value_attribute = "ALL"

            # what to do with those numbers. mean, median, ... default is "ALL", saves them all
            if(hasattr(self, "func") == False):
                self.func = "ALL"

                
            # Sets the interval where the values are recorded. 1 for each time step, 0 for never
            if(hasattr(self, "Dt") == False):
                self.Dt = 1
                
            # if end is true, it saves the state at the end of the simulation (might double up the last datapoint)
            if(hasattr(self, "end") == False):
                self.end = False
                
            # how do you want to call this item?
            if(hasattr(self, "label") == False):
                self.label = "Default_name_behavior"

            
        elif(self.TYPE == "DIS"):
            
            # name of the attribute,
            if(hasattr(self, "name_attribute") == False):
                self.name_attribute = "health_status"

            # value of the attribute, (for example in SIR, if you only want I set value_attribute = 1)    
            if(hasattr(self, "value_attribute") == False):
                self.value_attribute = "ALL"

            # what to do with those numbers. mean, median, ... default is "ALL", saves them all
            if(hasattr(self, "func") == False):
                self.func = "ALL"

            # Sets the interval where the values are recorded. 1 for each time step, 0 for never
            if(hasattr(self, "Dt") == False):
                self.Dt = 1
                
            # if end is true, it saves the state at the end of the simulation (might double up the last datapoint)
            if(hasattr(self, "end") == False):
                self.end = False
                
            # how do you want to call this item?
            if(hasattr(self, "label") == False):
                self.label = "Default_name_disease"
                
            if(hasattr(self, "tvec") == False):
                self.tvec = []

    def _print(self):
        for i in dir(self):
            if((i[0] != "_")):
                print(str(i) + " = " + str(getattr(self, i)))
    def _list_names(self):
        res = []
        for i in dir(self):
            if((i[0] != "_")):
                res.append(str(i))
        return res



def def_PS(target,file):
    
    """
    this creates the parameters for the Simulations.
    PS = def_PS(n=200, T=100, replicas=1, order = "DBD")
    n= number of node in the network
    T = number of simulation steps
    replicas = number of repetition of the simulation to take better averages
    order = update order. Default is Disease-Behavior-Disease.
            the options are: 
                D -- Disease^n
                DBD     --- (Disease Behavior)^n
                BDB     --- (Behavior Disease)^n
                BDD     --- precalculate Behavior (Disease)^n
                DBB     --- precalculate Disease (Behavior)^n
    """
    
    PS = import_par_json(target,file)
    if(hasattr(PS,"n") == False):
        #set default number of nodes:
        setattr(PS,"n",1000)
    if(hasattr(PS,"T") == False):
        #set default evolution steps:
        setattr(PS,"T",1000)
    if(hasattr(PS,"replicas") == False):
        #set default number of replicas per set of parameters:
        setattr(PS,"replicas",1)
    if(hasattr(PS,"order") == False):
        #set default order of simulation:
        setattr(PS,"order", "D")
    return(PS)


#def def_PD(network = "ER", link_m=6, name = "SIS", N_pat_zero = 1, 
#           S2I = 0, S2R = 0, S2V = 0, 
#           I2R = 0, I2S = 0, I2V = 0, 
#           R2S = 0, R2I = 0, R2V = 0, 
#           V2S = 0, V2I = 0, V2R = 0,
#           WS_p = 0):
def def_PD(target,file):
    
    PD = import_par_json(target,file)
    
    # adjust the network:
    
    
    if(hasattr(PD,"network") == False):
        #set default number of nodes:
        setattr(PD,"network","ER")
    if(PD.network == "ER"):
        if(hasattr(PD,"link_m") == False):
            #set default number of nodes:
            setattr(PD,"link_m",6)
            
            
    # Adjust the disease
    
    if(hasattr(PD,"N_pat_zero") == False):
        setattr(PD,"N_pat_zero",1)    

    if(hasattr(PD,"name") == False):
        #set default evolution steps:
        setattr(PD,"name","Sb0IS")
    
    if(PD.name  == "Sb0IS"):
        if(hasattr(PD,"S2I") == False):
            setattr(PD,"S2I",0.1)
        if(hasattr(PD,"S2I") == False):
            setattr(PD,"S2I",0.01)    
    
    IC = Parameters()
    if(PD.IC["location"] == "random"):
        IC.location = "random"
        IC.N_pat_zero = PD.IC["N_pat_zero"]
    elif(PD.IC["location"] == "exact"):
        IC.location = "exact"
        PD.IC.indx = PD.IC["indx_pat_zero"]
    elif(PD.IC["location"] == "central"):
        IC.location = "central"
        IC.indx = np.floor(PD.Lx*PD.Ly/2).astype(int)

    PD.IC = IC   # now it's not a dictionary anymore.
    
    
    
    return(PD)

    


def def_PB(target,file):
    
    PB = import_par_json(target,file)
    
    
    if(hasattr(PB,"network") == False):
        setattr(PB,"network","ER")

    if(PB.network == "ER"):
        if(hasattr(PB.network,"link_m") == False):
            setattr(PB.network,"link_m",10)
    
    if(hasattr(PB,"name") == False):
        setattr(PB,"name","No_B")
            
    if(PB.name == "Discrete"):
        if(hasattr(PB,"list_values") == False):
            setattr(PB,"list_values",[0,1])

            
    if(hasattr(PB, "CLUSTER")):
        CLUSTER = Parameters()
        for el in PB.CLUSTER:
            setattr(CLUSTER,el,PB.CLUSTER[el])
        PB.CLUSTER = CLUSTER

    
    if(hasattr(PB, "IC")):
        IC = Parameters()
        
        for i in PB.IC:
            
            temp = Parameters()

            
            if (PB.IC[i]["name"] == "discrete"):
                temp.name = "discrete"
                temp.values =  PB.IC[i]["values"]
                temp.weights = PB.IC[i]["weights"]
            
            elif(PB.IC[i]["name"] == "uniform"):
                temp.name = "uniform"
                temp.low = PB.IC[i]["low"]
                temp.up =  PB.IC[i]["up"]
            
            elif(PB.IC[i]["name"] == "beta"):
                temp.name = "beta"
                temp.beta_a = PB.IC[i]["beta_a"]
                temp.beta_b = PB.IC[i]["beta_b"]
            
            elif(PB.IC[i]["name"] == "constant"):
                temp.name = "constant"
                temp.value = PB.IC[i]["value"]            
                
            setattr(IC,i,temp)
        PB.IC = IC   # now it's not a dictionary anymore.
    
    
    return(PB)
    
    
    
def def_PR(target,file):
    # parameters for the recording
    PR = Parameters()
    PR.savings = import_par_json_savs(target, file)
    
    
    
    return(PR)




def import_par_json(FIRST_LABEL, NAMEFILE):
    with open(NAMEFILE) as f:
        DIC = json.load(f)

    if FIRST_LABEL in DIC:
        RES = Parameters()
        for el in DIC[FIRST_LABEL]:
            setattr(RES,el,DIC[FIRST_LABEL][el])
    else:
        RES = Parameters()
        print("WARNING!!! no " + FIRST_LABEL)
    return RES


def import_par_json_savs(FIRST_LABEL, NAMEFILE):
    
    with open(NAMEFILE) as f:
        DIC = json.load(f)
    #print(DIC)
    RES = []
    if FIRST_LABEL in DIC:
        for i in range(len(DIC[FIRST_LABEL])):
            # A = DIC[FIRST_LABEL][str(i+1)]
            # print(A)
            # print(type(A))
            T = SAVE(DIC[FIRST_LABEL][str(i+1)])
            #print(T)
            #print(T.end)
            T.end = bool(T.end)
            RES.append(T)
    else:
        RES = []
        print("WARNING!!! no SAVINGS")
    return RES





      

# the outcome of each replica will be saved here.

def init_result(Savings):
    
    result = Parameters()
    
    for key in Savings:
        
        name = key.label
        setattr(result, name, Parameters())
        setattr(getattr(result, name), "RES", [])
        setattr(getattr(result, name), "Tvec", [])
        setattr(getattr(result, name), "PAR", key)   
    
    return result




    
    
def iscondition(PR,condition,LIST_save):
    for val in LIST_save:
        if(val.func ==condition):
            return True
    return False



#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import igraph as ig

#from numba import jit

def grid_layout(rows, cols):
    layout = []
    for row in range(rows):
        for col in range(cols):
            layout.append((col, rows - row - 1))
    return layout


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import igraph as ig
import random as rn
import pathlib
#from numba import jit

from Parameter_Importer import *



def my_filtering_target(g, node, xyz_neighbors, vertex_attribute_name, vertex_value_required):
    
    res = []
    

    for v in g.vs[node[xyz_neighbors]]:
        if(v[vertex_attribute_name] == vertex_value_required):
            res.append(v)
    return res

        
def normalize(g,PS,PR):
    PR.RES=PR.RES/PS.n/PS.replicas
    PR.RES_B=PR.RES_B/PS.replicas
    
def calc_polarization_var(o):
    o = np.array(o)
    #print(o)
    M = np.ones([len(o), len(o)]) * o
    result = np.std(np.abs(M*o - (M*o).T) )**2
    
    # o = np.array(o)
    # diff = np.array([np.abs(o[j] - o[i]) for i in range(len(o)) for j in range(i+1, len(o))])
    
    return result #np.var(diff)

def calc_polarization_avg(o):
    o = np.array(o)
    #print(o)
    M = np.ones([len(o), len(o)]) * o
    result = np.average(np.abs(M*o - (M*o).T) )
    
    # o = np.array(o)
    # diff = np.array([np.abs(o[j] - o[i]) for i in range(len(o)) for j in range(i+1, len(o))])
    
    return result #np.var(diff)





#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import igraph as ig
import random as rn
from scipy.stats import beta
import pathlib
#from numba import jit
from Save_Functions import *


from Parameter_Importer import *
from Clustering import *


def generate_adj_matrix(PS,PX, graph_to_copy = -1):
    
    if (PX.network == "ER"):
        er_p=min(1,PX.link_m/(PS.n))

        g_xyz = ig.Graph.Erdos_Renyi(n=PS.n, p=er_p)
    elif (PX.network == "Lattice"):
        Lx = PX.Lx
        Ly = PX.Ly
        if(Lx*Ly != PS.n):
            input("ERROR! the size of the lattice graph is not the one you wanted")
            PS.n=Lx * Ly
        g_xyz=ig.Graph.Lattice(dim=[Lx, Ly], circular=False)
    elif (PX.network == "WS"):
        #g_xyz = ig.Graph.Watts_Strogatz(ndim=PX.ndim_WS, size=PX.size_WS, nei=PX.nei_WS, p=PX.WS)
        g_xyz = ig.Graph.Watts_Strogatz(dim=1, size=PS.n, nei=2, p=PX.WS_p)
    elif (PX.network == "same"):
        if(graph_to_copy == -1):
            input("NO COPY NET GIVEN")
        else:
            g_xyz = graph_to_copy.copy()
    else:
        print("GRAPH NOT IMPLEMENTED YET! NUUUUUU")

    return g_xyz    
        
def initialize_disease(g,PS,PD):
        
    # extract N patient zero from where the disease starts: 
    # and set their status at 2, while all the other nodes are set at 1.
    #print(PD.name)

    if((PD.name == "SIRB")|(PD.name == "SIR_B_recovery")|(PD.name == "No_D")):
        
        if(PD.IC.location == "random"):
            indx = (rn.sample(range(0, PS.n), PD.IC.N_pat_zero))
            g.vs["health_status"]=1
            g.vs[indx]["health_status"]=2
        elif(PD.IC.location == "center"):
            g.vs["health_status"]=1
            indx = np.floor(31*31/2) # this is the center of the graph, it makes sense only for the lattice
            g.vs[indx]["health_status"]=2
        else:
            g.vs["health_status"]=1
            g.vs[PD.IC.indx]["health_status"]=2

    else:
        print("WARNING! disease unknown. standard initalization")
        indx = np.array(rn.sample(range(0, PS.n), PD.N_pat_zero))
        g.vs["health_status"]=1
        g.vs[indx]["health_status"]=2
        
        
        
        
        
        
        
        
        
        

def assign_IC(IC,target,N,
              d_name = "dummy", d_value = 0, d_values = [0,1], 
              d_weights = [0.5, 0.5], d_beta_a = 0.1, 
              d_beta_b = 0.1, d_low=0, d_high = 1):

        #print("target =" + target)
    
        if(hasattr(IC,target)):
            A = getattr(IC,target)
            #print("A.name = " + A.name)
            if(A.name == "uniform"):
                #print("uniform" + str(A.low) + " --- " + str(A.high))
                return np.random.uniform(low = A.low, high = A.high, size = N)
            if(A.name == "discrete"):
                #print("discrete" + str(A.values) + " --- " + str(A.weights))
                return np.random.choice(A.values, N, p = A.weights)
            
            if(A.name == "constant"):
                #print("constant" + str(A.value))
                return (np.ones(N)*A.value)
            if(A.name == "beta"):
                #print("beta" + str(A.beta_a) + " --- " + str(A.beta_b))
                return np.array(beta.rvs(A.beta_a, A.beta_b, size=N))
        else:
            #print("default")
            if(d_name == "uniform"):
                return np.random.uniform(low = d_low, high = d_high, size = N)
            if(d_name == "discrete"):
                return np.random.choice(d_values, N, p = d_weights)
            if(d_name == "constant"):
                return (np.ones(N)*d_value)
            if(d_name == "beta"):
                return np.array(beta.rvs(d_beta_a, d_beta_b, size=N))
            
        
def initialize_behavior(g,PS,PB,CONSTANTS):
    
    # Purely behavioral models

    if(PB.name == "Average"):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n, d_name="uniform", d_low = 0, d_high = 1) # def = uniform [0,1]
            
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)     #initialize the initial behavior
    
    elif(PB.name == "Reinforcement"):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n, d_name="uniform", d_low = 0, d_high = 1) # def = uniform [0,1]
            
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)     
    
    elif(PB.name == "Repulsive"):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n, d_name="uniform", d_low = 0, d_high = 1) # def = uniform [0,1]
            
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)     
    elif(PB.name == "BCM"):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n, d_name="uniform", d_low = 0, d_high = 1) # def = uniform [0,1]            
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)     
    elif(PB.name == "Negative_Influence"):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n, d_name="uniform", d_low = 0, d_high = 1) # def = uniform [0,1]            
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)     
    





    elif((PB.name == "Xia-like-Heaviside")):
        
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S", PS.n, d_name="uniform", d_low = 0, d_high = 1) # def = uniform [0,1] 
            CONSTANTS.P_P = assign_IC(PB.IC,"P_P", PS.n, d_name="uniform", d_low = 0, d_high = 1) # def = uniform [0,1]
            
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)     
        g.vs['peer_pressure']   = np.random.permutation(CONSTANTS.P_P)     
        
        
    elif((PB.name == "Xu2022")|(PB.name == "Xu2022_light")):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n, d_name="constant", d_value = 0) # def = 0 
            
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)     
        
        
        
    elif(PB.name == "PIRES"):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n, d_name="constant", d_value = 0) # standard here constant = 0
            CONSTANTS.P_P = assign_IC(PB.IC,"P_P",PS.n, d_name="uniform", d_low = 0, d_high = 1) # def = uniform [0,1]
            
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)     #initialize the initial behavior
        g.vs['peer_pressure']   = np.random.permutation(CONSTANTS.P_P)     #initialize the initial behavior


    
    elif(PB.name == "XIA_2013"):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n, d_name = "discrete", d_values=[0,1], d_weights = [0.5, 0.5]) # def = discrete [0,1]
        
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)
        g.vs["new_state_B"] = g.vs['behavior_status']
        g.vs["health_status"] =g.vs['behavior_status']
        g.vs["L_avg_D"] = np.ones(PS.n)
    
    
    elif(PB.name == "No_B"):
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = assign_IC(PB.IC,"B_S",PS.n)
            #print(CONSTANTS.B_S)
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)
        
    else:
        print("WARNING! behaviour unknown. standard initalization, fingers crossed")
        print(PB.name)
        if(CONSTANTS.flag==0):
            CONSTANTS.flag = 1
            CONSTANTS.B_S = np.random.choice(PB.list_values, size=PS.N)
        
        g.vs['behavior_status'] = np.random.permutation(CONSTANTS.B_S)


    
    
def init_graph(PS,PB,PD,PR):
    # Generate the disease-graph:
    g_dis= generate_adj_matrix(PS,PD)
    g_beh= generate_adj_matrix(PS,PB,g_dis)
    # Create an attribute for each link type. 
    # To distinguish those who are disease connected from those who are behavior connected
    g_dis.es["dis_linked"]=True
    g_beh.es["beh_linked"]=True    
    # label the indices of the two networks
    ca = 0
    for v in g_dis.vs:
        v['index_dis'] = ca
        ca = ca+1
        
    cb = 0
    for v in g_beh.vs:
        v['index_beh'] = cb
        cb = cb+1
    
    
    #precalculate the index of the neighbors of each node. And the number of neighbors.
    for v in g_dis.vs:
        neigh=v.neighbors()
        v['dis_neighbors']=[v2['index_dis'] for v2 in neigh]
        v['L_dis_neighbors'] = len(v['dis_neighbors'])
    
    for v in g_beh.vs:
        neigh=v.neighbors()
        v['beh_neighbors']=[v2['index_beh'] for v2 in neigh]
        v['L_beh_neighbors'] = len(v['beh_neighbors'])
    
    
    g=g_dis.union(g_beh) # this maybe too ^^^^^
    
    # Assign indexes to the merged graph  # this maybe too ^^^^^^^
    count = 0
    for v in g.vs:
        v['indx'] = count
        count = count+1
    
    
    return g
        

def init_sim(g, CONSTANTS, PS,PB,PD,PR,result):
    
    
    initialize_disease(g,PS,PD)
    initialize_behavior(g,PS,PB,CONSTANTS)
    
    
    
    if(hasattr(PB, "CLUSTER")):
        Hom_history_evo = cluster_behavior(g,PB.CLUSTER,debug = True)        
        if(hasattr(result,"homophily_curve") == True):
            result.homophily_curve.RES = Hom_history_evo        

    # Merges the two networks