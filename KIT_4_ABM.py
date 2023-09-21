#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import igraph as ig
import random as rn
import pathlib
from Clustering import *
from Init_and_Dynamics import *
from Save_Functions import *
from Utilities import *
import json

def import_parameters(namefile):
    # namefile is a txt file in the JSON format.
    # first I import everything into data, then I split data into P_lay, P_dyn, P_sim, P_rec

    with open(namefile) as f:
        data = json.load(f)

    P_lay = data["P_layer"]
    P_dyn = data["P_dynamic"]
    P_sim = data["P_simulations"]
    P_rec = data["P_recordings"]

    return P_lay, P_dyn, P_sim, P_rec


def run_sim(P_layer, P_dynamic, P_simulations, P_recordings):

    # initialize the graph creating all the needed layers. for each layer create the network
    G = init_graph(P_layer) # G is a list of graphs.

    # imprints the rules on the graph and returns a list of rules, and of layers to execute for the corresponding rule.
    list_of_rules = imprint_rules(G, P_dynamic)
    
    #print("The graph is now completed. now running the simulations...")
    run_temporal_evolution(G, list_of_rules, P_simulations, P_recordings)

#  ██ ███    ██ ██ ████████                     ███    ██ ███████ ████████ 
#  ██ ████   ██ ██    ██                        ████   ██ ██         ██    
#  ██ ██ ██  ██ ██    ██        █████ █████     ██ ██  ██ █████      ██    
#  ██ ██  ██ ██ ██    ██                        ██  ██ ██ ██         ██    
#  ██ ██   ████ ██    ██                        ██   ████ ███████    ██    

def init_graph(P_lay):
    G = []
    for i in range(P_lay["N_layers"]):
        p_lay_i = P_lay["Layer_" + str(i)]
        N = P_lay["N_nodes"]

        if(p_lay_i["type"] == "ER_m"):
            g = ig.Graph.Erdos_Renyi(n=N, p=p_lay_i["link_m"]/(N-1))
        elif(p_lay_i["type"] == "ER_p"):
            g = ig.Graph.Erdos_Renyi(n=N, p=p_lay_i["p"])
        elif(p_lay_i["type"] == "file"):
            #if p_lay_i has not the field "format" try to load the file with igraph.Load without forcing the format:
            if("format" not in p_lay_i):
                g = ig.Graph.Load(p_lay_i["filename"])
            else:
                g = ig.Graph.Read_Ncol(p_lay_i["filename"], format = p_lay_i["format"])
        elif(p_lay_i["type"] == "Lattice"):
            if(p_lay_i["Lx"]*p_lay_i["Ly"] != N):
                print("GRAPH SIZE WARNING: Lx*Ly != N: " + str(p_lay_i["Lx"]*p_lay_i["Ly"]) + " != " + str(N)) 
            g = ig.Graph.Lattice(dim=[p_lay_i["Lx"], p_lay_i["Ly"]], circular=False)
        elif(p_lay_i["type"] == "WS"):
            if(np.power(p_lay_i["L"],p_lay_i["D"]) != N):
                print("GRAPH SIZE WARNING: L^D != N: " + str(np.power(p_lay_i["L"],p_lay_i["D"])) + " != " + str(N))
            g = ig.Graph.Watts_Strogatz(dim=p_lay_i["D"], size=p_lay_i["L"], nei=p_lay_i["NFN"], p=p_lay_i["P"])
        else:
            print("GRAPH: " + p_lay_i["NETWORK"]["type"] + " NOT IMPLEMENTED YET! NUUUUUU")

        G.append(g)
    return G


#  ██ ███    ███ ██████  ██████  ██ ███    ██ ████████     ██████  ██    ██ ██      ███████ ███████ 
#  ██ ████  ████ ██   ██ ██   ██ ██ ████   ██    ██        ██   ██ ██    ██ ██      ██      ██      
#  ██ ██ ████ ██ ██████  ██████  ██ ██ ██  ██    ██        ██████  ██    ██ ██      █████   ███████ 
#  ██ ██  ██  ██ ██      ██   ██ ██ ██  ██ ██    ██        ██   ██ ██    ██ ██      ██           ██ 
#  ██ ██      ██ ██      ██   ██ ██ ██   ████    ██        ██   ██  ██████  ███████ ███████ ███████ 
#                                                                                                   

def imprint_rules(G,P_dyn):

    LotR = [] # list of the rules
    for i in range(P_dyn["N"]):
        # for each dynamic i, read the parameters and initialize the rule
        P_dyn_i = P_dyn["Dynamic_" + str(i)]

        
        single_rule = init_single_rule[P_dyn_i["func"]](P_dyn_i, G)
        
        
        LotR.append(single_rule)

    return LotR # list of rules




def run_temporal_evolution(G, list_of_rules, P_simulations, P_recordings):
    # G is the list of graphs-layers (for now a single graph)
    # list_of_rules is the list of rules
    # P_recordings is the dictionary with the parameters for the recordings


    # initialize the recordings -- to be updated, since I now want to save everything as h5 files
    # L_REC_0 is the list of recordings to be done BEFORE the simulations begin
    # L_REC is the list of recordings to be done DURING the simulations
    # L_REC_1 is the list of recordings to be done AFTER the simulations end

    L_REC_0, L_REC, L_REC_1 = init_recordings(P_recordings, P_simulations["T"])

    # save the state of the system BEFORE the simulations begin
    for record0 in L_REC_0:
        single_save(G, record0, -2)     # -1 means that the time step is before the simulations begin

    for internal_tick in range(P_simulations["T"]):
        #print ("Time step: " + str(i) + "/" + str(P_simulations["T"]))

        # for each time step internal_tick, run the simulation
        for rule in list_of_rules:
            single_update(G, rule)
        
        # save the state of the system DURING the simulations (if the record delay is a multiple of the internal_tick, allows to save the state every N steps)
        for record in L_REC:
            if(internal_tick%record["DT"] == 0):
                single_save(G, record, internal_tick, Delta_T = record["DT"])
    
    # save the state of the system AFTER the simulations end
    for record1 in L_REC_1:
        single_save(G, record1, -1)     # -1 means that the time step is after the simulations end





def single_update(G, rule):
    rule_name = rule["rule"]
    try:
        function = update_dictionary[rule_name]
        function(G, rule)  # Call the function
    except KeyError:
        print("I cannot update the following rule:" + rule_name)
        
        print(type(rule_name))
        print("ERROR: single_update of rule(" + rule + ") NOT IMPLEMENTED YET! NUUUUUU")
    pass


# %%
