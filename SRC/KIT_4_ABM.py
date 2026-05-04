#!/usr/bin/env python
# coding: utf-8

#import modules
import numpy as np
import igraph as ig
import sys
from copy import deepcopy
from difflib import get_close_matches    #for suggesting update_function_keys
import json

#import own modules
from Load_Models import *   #this script runs and imports all functions for Health-Behavior models defined in the directory models
from Save_Functions import *   #functions for extracting information from the simulation
from Write_functions import *    #functions to save write extracted information to disk
from utilities.network import *


class Global_Vars:
    def __init__(self):
        self.stop_condition = False

#  █████      ██     █████      ██     ██   ██  ██████   ██████   ██████   █████     ████   
#  ██  ██    ████    ██  ██    ████    ███ ███  ██         ██     ██       ██  ██   ██  ██  
#  ██  ██   ██  ██   ██  ██   ██  ██   ███████  ██         ██     ██       ██  ██   ██      
#  █████    ██████   █████    ██████   ██ █ ██  ████       ██     ████     █████     ████   
#  ██       ██  ██   ████     ██  ██   ██   ██  ██         ██     ██       ████         ██  
#  ██       ██  ██   ██ ██    ██  ██   ██   ██  ██         ██     ██       ██ ██    ██  ██  
#  ██       ██  ██   ██  ██   ██  ██   ██   ██  ██████     ██     ██████   ██  ██    ████   

def import_parameters(namefile):
    # namefile is a json file in the JSON format.
    # first I import everything into data, then I split data into P_net, P_dyn, P_sim, P_rec

    with open(namefile) as f:
        data = json.load(f)

    P_net = data.get("P_network",data.get("P_layer"))
    P_dyn = data["P_dynamic"]
    P_sim = data["P_simulations"]
    P_rec = data.get("P_record",data.get("P_recordings"))

    return P_net, P_dyn, P_sim, P_rec

def reset_param(P_net_original,P_sim_original,P_dyn_original,P_rec_original):
    P_net = deepcopy(P_net_original)
    P_dyn = deepcopy(P_dyn_original)
    P_sim = deepcopy(P_sim_original)
    P_rec = deepcopy(P_rec_original)
    return P_net, P_sim, P_dyn, P_rec

  
#  ██  ██   █████    ████       ██     ██████   ██████  
#  ██  ██   ██  ██   ██ ██     ████      ██     ██      
#  ██  ██   ██  ██   ██  ██   ██  ██     ██     ██      
#  ██  ██   █████    ██  ██   ██████     ██     ████    
#  ██  ██   ██       ██  ██   ██  ██     ██     ██      
#  ██  ██   ██       ██ ██    ██  ██     ██     ██      
#   ████    ██       ████     ██  ██     ██     ██████  
#                                                                                                                                                    

def init_update_rules(G,P_dyn,global_var):

    LotR = [] # list of the rules
    for i in range(P_dyn["N"]):
        # for each dynamic i, read the parameters and initialize the rule
        P_dyn_i = P_dyn["Dynamic_" + str(i)]

        try:
            init_function = init_fct_dict[P_dyn_i["func"]]
        except KeyError:
            print(f"Invalid init function key: {P_dyn_i["func"]}")
            suggestion = get_close_matches(P_dyn_i["func"], init_fct_dict.keys(), n=1)
            if suggestion:
                print(f"Did you mean {suggestion[0]} ?")
            exit()

        P_rule = init_function(P_dyn_i, G, global_var)
        #here the init function is called    (referenced in init_fct_dict with the key P_dyn_i["func"])
        #That function initializes the dynamics
        # and returns a    P_rule    , a dictionary with all parameters necessary for the update function

        LotR.append(P_rule)

    return LotR # list of rules


def single_update(G, P_rule, global_var = Global_Vars()):
    update_fct_name = P_rule["func"]

    try:
        update_fct = update_fct_dict[update_fct_name]
    except KeyError:
        print(f"Invalid update function key: {update_fct_name}")
        suggestion = get_close_matches(update_fct_name, update_fct_dict.keys(), n=1)
        if suggestion:
            print(f"Did you mean {suggestion[0]} ?")

    update_fct(G, P_rule, global_var)  # Call the function


#  █████    ██  ██   ██  ██             ████     ████    ██   ██ 
#  ██  ██   ██  ██   ███ ██            ██  ██     ██     ███ ███ 
#  ██  ██   ██  ██   ██████            ██         ██     ███████ 
#  █████    ██  ██   ██████             ████      ██     ██ █ ██ 
#  ████     ██  ██   ██ ███                ██     ██     ██   ██ 
#  ██ ██    ██  ██   ██  ██            ██  ██     ██     ██   ██ 
#  ██  ██    ████    ██  ██             ████     ████    ██   ██ 


def run_sim(P_network, P_dynamic, P_simulations, P_record, initialized_recording = None, return_G = False):
    #initialize the graph, read in all simulation rules, simulate and return data

    # initialize the graph creating all the needed layers. for each layer create the network
    Graphs = init_graph(P_network) # G is a list of graphs.

    #initialize the container object global_var
    global_var = Global_Vars()
    #global_var contains variables that should be available everywhere
    global_var.functions = Global_Vars()
    #global_var.functions contains fucntions that should be available everywhere

    # initialize the dynamic on the graph and return a list of rules for the updating function.
    update_rules = init_update_rules(Graphs, P_dynamic,global_var)

    #initialize the results object and if necessary the recording lists
    if initialized_recording == None:
        L_REC_begin, L_REC_during, L_REC_after = init_recording(P_record, P_simulations["T"], P_network)
    else:
        #initialized_recording allows to initialize everything outside of simulation - and skip doing it all over again here
        #for saving to h5 you need to know the shape of everything beforehand, and init_recording finds out that information
        # also useful when doing multiple trials or parameter sweeps:  all L_REC lists stay the same
        L_REC_begin, L_REC_during, L_REC_after = initialized_recording

    results = init_results(P_record)
    #results need to be initialized for every trial so the old ones are deleted

    global_var.current_timestep = 0
    # save the state of the system BEFORE the simulations begins
    for P_rec_i in L_REC_begin:
        #P_rec_i contains all parameters for saving one single variable before the simulation
        single_save(Graphs, P_rec_i, results, global_var, internal_tick = 0)     # 0 means that the time step is before the simulations begin
        #writes data to              results

    #simulate
    for internal_tick in range(1,P_simulations["T"]+1):
        global_var.current_timestep = internal_tick

        if not global_var.stop_condition:
            # for each time step, advance the simulation one increment
            for P_rule in update_rules:
                single_update(Graphs, P_rule, global_var)
        
            # save the state of the system DURING the simulation
            # with the possibility to save the state every N steps via P_rec_i["DT"]
            for P_rec_i in L_REC_during:
                if(internal_tick%P_rec_i["DT"] == 0):
                    single_save(Graphs, P_rec_i, results, global_var, internal_tick = internal_tick)
                    #writes data to              results

        if global_var.stop_condition:
            for P_rec_i in L_REC_during:
                batch_save(Graphs, P_rec_i,results, global_var, internal_tick = internal_tick, T = P_simulations["T"])
                #writes lots of data to    results
            break

    # save the state of the system AFTER the simulations end
    for P_rec_i in L_REC_after:
        single_save(Graphs, P_rec_i, results, global_var, internal_tick = -1)     # -1 means that the time step is after the simulations end
        #writes data to              results

    convert_results_to_float(results)    

    if return_G:
        return results, Graphs
    return results









