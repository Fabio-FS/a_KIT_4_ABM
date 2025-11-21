#!/usr/bin/env python
# coding: utf-8


import numpy as np
import igraph as ig
from Init_and_Dynamics import *
from Save_Functions import *
from Write_functions_csv import *
import sys
from copy import deepcopy

import json

class Global_Vars:
    def __init__(self):
        self.stop_condition = False


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

def run_sim(P_network, P_dynamic, P_simulations, P_record, return_G = False):

    # initialize the graph creating all the needed layers. for each layer create the network
    Graphs = init_graph(P_network) # G is a list of graphs.

    global_var = Global_Vars()
    #global_var contains variables that all functions should have access to in principle
    #for the moment it's only stop_condition

    # initializes the dynamic on the graph and returns a list of rules for the updating function.
    list_of_rules = init_rules(Graphs, P_dynamic,global_var)
    
    Data = simulate_and_return_data(Graphs, list_of_rules, P_simulations, P_record,global_var)

    if return_G:
        return Data, Graphs
    return Data

# ██  ██   ██████   ██████   ██   ██   ████    █████    ██  ██  
# ███ ██   ██         ██     ██   ██  ██  ██   ██  ██   ██ ██   
# ██████   ██         ██     ██   ██  ██  ██   ██  ██   ████    
# ██████   ████       ██     ██ █ ██  ██  ██   █████    ███     
# ██ ███   ██         ██     ███████  ██  ██   ████     ████    
# ██  ██   ██         ██     ███ ███  ██  ██   ██ ██    ██ ██   
# ██  ██   ██████     ██     ██   ██   ████    ██  ██   ██  ██  


def init_graph(P_net):
    G = []
    for i in range(P_net["N_layers"]):
        P_layer = P_net["Layer_" + str(i)]
        N = P_net["N_nodes"]

        if(P_layer["type"] == "ER_m"):
            g = ig.Graph.Erdos_Renyi(n=N, p=P_layer["link_m"]/(N-1))
        elif(P_layer["type"] == "ER_p"):
            g = ig.Graph.Erdos_Renyi(n=N, p=P_layer["p"])
        elif(P_layer["type"] == "file"):
            #if P_layer has not the field "format" try to load the file with igraph.Load without forcing the format:
            if("format" not in P_layer):
                g = ig.Graph.Load(P_layer["filename"])
            else:
                g = ig.Graph.Read_Ncol(P_layer["filename"], format = P_layer["format"])
        elif(P_layer["type"] == "Lattice"):
            if(P_layer["Lx"]*P_layer["Ly"] != N):
                print("GRAPH SIZE WARNING: Lx*Ly != N: " + str(P_layer["Lx"]*P_layer["Ly"]) + " != " + str(N)) 
            g = ig.Graph.Lattice(dim=[P_layer["Lx"], P_layer["Ly"]], circular=P_layer.get("circular",False), nei = P_layer["nei"])
        elif(P_layer["type"] == "Moore_Lattice"):
            #a lattice with Moore neighborhood, i.e. 8 neighbors
            #a Cellular automaton
            #number of agents
            N = P_layer["Lx"]*P_layer["Ly"]
            #set up adjacency matrix
            A = np.zeros(shape = (N,N))

            #horizontal: +1
            #vertical:   +x  (conflict with horiz. if x=1)
            #NW-SE:      +x+1
            #NE-SW:      +x-1 (conflict with horiz. if x=2)

            if P_layer["Lx"] == 1:
                #special case: line graph
                A += np.diag(np.ones(shape = N-1), 1 )    #→→→→ vertical = horizontal conections +1 (in adj matrix: one off the diagonal)
                if P_layer["circular"]:
                    A[0,N-1] = 1
                A += A.T
                 
            elif P_layer["Lx"] == 2:
                #special treatment because: ↙↙ NE-SW connections are same distance as →→ horizontal connections
                for i in range(N):
                    if i == 0:   #first row
                        A[i,:4] = np.array([0,1,1,1])
                    elif i==1:    #second row
                        A[i,:4] = np.array([1,0,1,1])
                    elif i==N-2:    #second to last row
                        A[i,-4:] = np.array([1,1,0,1])
                    elif i==N-1:    #last row
                        A[i,-4:] = np.array([1,1,1,0])
                    elif i%2 == 0:                             #every second row
                        print(i,2*(i//2)-2)
                        A[i,2*(i//2)-2:2*(i//2)-2+6] = np.array([1,1,0,1,1,1])
                    else:                                      #every other second row
                        A[i,2*(i//2)-2:2*(i//2)-2+6] = np.array([1,1,1,0,1,1])
                if P_layer["circular"] == True:
                    A[N-1,1] = 1    #connnecting last agent with 2nd
                    A[N-2,0] = 1    #connecting 2nd to last agent with first
                    A[N-1,0] = 1    #connnecting last agent with 1st
                    A[N-2,1] = 1    #connecting 2nd to last agent with 2nd
                A += A.T
                A = np.clip(A,a_min=0,a_max=1)

            else:   #xdim is > 2
                Lx = P_layer["Lx"]
                #examples always N=9, Lx=3

                #adding 1s to adjacency matrix
                #including few ones that have to be deleted again
                A += np.diag(np.ones(shape = N-1), 1 )                            #→→→→ horizontal conections +1 (in adj matrix: one off the diagonal)
                A[np.arange(0, N-Lx,   1)   ,   np.arange(Lx,   N, 1)]   =   1    #↓↓↓↓  vertical    +x
                #e.g.    0,1,2,3,4,5                  3,4,5,6,7,8      
                A[np.arange(0, N-Lx-1, 1)   ,   np.arange(Lx+1, N, 1)]   =   1    # ↘↘↘↘  NW-SE       +x+1
                #e.g.    0,1,2,3,4                    4,5,6,7,8
                A[np.arange(0, N-Lx+1, 1)   ,   np.arange(Lx-1,   N, 1)]   =   1    # ↙↙↙↙ NE-SW       +x-1
                #e.g.    0,1,2,3,4,5,6                  2,3,4,5,6,7,8

                #making adjacency matrix symmetric
                A += A.T

                #corrections:
                #1. →→→→ horizontal correction, removing "connections through line break"
                A[ np.arange(Lx-1, N-1, Lx)   ,   np.arange(Lx, N, Lx)]  = 0      
                #e.g.        2,5                              3,6
                #
                #2. ↘↘↘↘ NW-SE correction;     removing connection from rightmost in one line to leftmost 2 lines below
                A[ np.arange(Lx-1, N-Lx-1, Lx )   ,   np.arange(2*Lx,N,Lx )] = 0
                #e.g.        2                                  6
                #
                #3. ↙↙↙↙ NE-SW   correction:    removing connection between left-most and right-most in same line   
                if not P_layer["circular"]:
                    A[np.arange(0,N,Lx)   ,   np.arange(Lx-1,N+Lx-1,Lx)] = 0        
                    #e.g.     0,3,6                      2,5,8  
                    
                #removing these connections also in reverse direction
                A *= A.T

                #adding connections for circular lattices
                if P_layer["circular"] == True:

                    #1. ↓↓↓↓ adding vertical connections between first and last line
                    A[np.arange(0,Lx,1)     ,    np.arange(N-Lx,N,1)]    =   1
                    #e.g.       0,1,2                      6,7,8
                    #
                    #2. ↘↘↘↘ adding NW-SE connections on the right edge
                    A[np.arange(Lx-1,N,Lx)    ,   np.mod(np.arange(Lx,N+1,Lx),N)]   =   1
                    #e.g.       2,5,8                            3,6,9mod9=0
                    #
                    #3. ↙↙↙↙ adding NE-SW connections on right edge
                    A[np.arange(Lx-1,N,Lx)    ,   np.mod(np.arange(-Lx,N-Lx,Lx),N)]   =   1
                    #e.g.        2,5,8                           6,0,3   (-3,0,3)
                    #
                    #4. ↘↘↘↘ adding NW-SE connections in bottom row, except bottom-right
                    A[np.arange(N-Lx,N-1,1)    ,    np.arange(1,Lx,1)]   =  1
                    #e.g.        6,7                      1,2
                    #5. ↙↙↙↙ adding NW-SE connecs in bottom row, except bottom-left
                    A[np.arange(N-Lx+1,N,1)    ,    np.arange(0,Lx-1,1)]   =  1
                    #e.g.        7,8                           0,1

                    A += A.T
                    A = np.clip(A,a_min=0,a_max=1)
            A = np.asarray(A, dtype=int)
            g = ig.Graph.Adjacency(A, mode="undirected")
        elif(P_layer["type"] == "WS"):
            if(np.power(P_layer["L"],P_layer["D"]) != N):
                print("GRAPH SIZE WARNING: L^D != N: " + str(np.power(P_layer["L"],P_layer["D"])) + " != " + str(N))
            g = ig.Graph.Watts_Strogatz(dim=P_layer["D"], size=P_layer["L"], nei=P_layer["NFN"], p=P_layer["P"])
        elif(P_layer["type"] == "kRRG"):
            g = ig.Graph.K_Regular(n = N, k = P_layer["k"])
        elif(P_layer["type"] == "BA"):
            g = ig.Graph.Barabasi(n = N, m =  P_layer["m"])
        elif(P_layer["type"] == "2islands"):

            p = P_layer["p"]
            k = P_layer["k"]
            #island1 = ig.Graph.Erdos_Renyi(n=N//2, p=p_lay_i["k"]/(N//2-1))
            #island2 = ig.Graph.Erdos_Renyi(n=N-N//2, p=p_lay_i["k"]/((N-N//2)-1))

            membership = np.empty(N)
            membership[:N//2] = np.full(N//2,fill_value=0)
            membership[N//2:] = np.full(N-N//2,fill_value=1)

            gs0 = N//2
            gs1 = N - N//2
            success = False
            while success == False:

                #print(2*p*k/(N-1), 2*(1-p)*k/(N-1), p*k/(N-1) + (1-p)*k/(N-1))
                #print(np.log(N)/N)
                adj_matrix00 = np.random.choice([0,1],size = (gs0,gs0), p=[1-2*(1-p)*k/(N-1), 2*(1-p)*k/(N-1)])
                adj_matrix11 = np.random.choice([0,1],size = (gs1,gs1), p=[1-2*(1-p)*k/(N-1), 2*(1-p)*k/(N-1)])
                adj_matrix01 = np.random.choice([0,1],size = (gs0,gs1), p=[1-2*p*k/(N-1), 2*p*k/(N-1)])

                adj_matrix00 = np.tril(adj_matrix00, -1) + np.tril(adj_matrix00, -1).T
                adj_matrix11 = np.tril(adj_matrix11, -1) + np.tril(adj_matrix11, -1).T
                #adj_matrix01 = np.tril(adj_matrix01, -1) + np.tril(adj_matrix01, -1).T
                #adj_matrix10 = np.tril(adj_matrix10, -1) + np.tril(adj_matrix10, -1).T

                #print(adj_matrix00)

                adj_matrix = np.zeros((N,N))

                adj_matrix[:gs0,:gs0] = adj_matrix00
                adj_matrix[gs0:,gs0:] = adj_matrix11
                adj_matrix[gs0:,:gs0] = adj_matrix01
                #adj_matrix[:gs0,gs0:] = adj_matrix10
                adj_matrix[:gs0,gs0:] = adj_matrix01.T

                g =  ig.Graph.Adjacency(adj_matrix.tolist(), mode=ig.ADJ_UNDIRECTED)


                if not any(x == 0 for x in g.degree()):
                    success = True
            #print("sum of adj matrix",np.sum(adj_matrix))
            #print("column sum of adj matrix",np.sum(adj_matrix,axis=0))
            #print("row sum of adj matrix",np.sum(adj_matrix,axis=1))
            #vertex_color = [(0,0,0) for i in range(N)]
            #for i in range(N):
            #    if membership[i] == 0:
            #        vertex_color[i] = (1,0,0)
            #    else:
            #        vertex_color[i] = (0,1,0)
            #ig.plot(g,target="sample_network.png",vertex_color=vertex_color)
            g.vs["membership"] = membership
            #print(g.vs["membership"])
        else:
            print("GRAPH: " + P_layer["type"] + " not implemented yet")

        G.append(g)
    return G

#  ██ ███    ███ ██████  ██████  ██ ███    ██ ████████     ██████  ██    ██ ██      ███████ ███████ 
#  ██ ████  ████ ██   ██ ██   ██ ██ ████   ██    ██        ██   ██ ██    ██ ██      ██      ██      
#  ██ ██ ████ ██ ██████  ██████  ██ ██ ██  ██    ██        ██████  ██    ██ ██      █████   ███████ 
#  ██ ██  ██  ██ ██      ██   ██ ██ ██  ██ ██    ██        ██   ██ ██    ██ ██      ██           ██ 
#  ██ ██      ██ ██      ██   ██ ██ ██   ████    ██        ██   ██  ██████  ███████ ███████ ███████ 
#                                                                                                   

def init_rules(G,P_dyn,global_var):

    LotR = [] # list of the rules
    for i in range(P_dyn["N"]):
        # for each dynamic i, read the parameters and initialize the rule
        P_dyn_i = P_dyn["Dynamic_" + str(i)]
        P_rule = init_fct_dict[P_dyn_i["func"]](P_dyn_i, G, global_var)
        #here the init function is called (referenced in init_fct_dict with the key P_dyn_i["func"]
        #That function initializes the dynamics
        # and returns a P_rule, a dictionary with all parameters necessary for the update function

        LotR.append(P_rule)

    return LotR # list of rules

def simulate_and_return_data(Graphs, list_of_rules, P_simulations, P_record,global_var):
    # G is the list of graph-layers, each item is one graph
    # list_of_rules is a list of dictionaries, each containing parameters for one updating function
    # P_record is the dictionary with the parameters for recording

    np.set_printoptions(threshold=sys.maxsize)
    #necessary for networks with >999 agents.
    #otherwise numpy prints every array as [y_0, y_1, ..., y_n]

    L_REC_0, L_REC, L_REC_1, results = init_recording(P_record, P_simulations["T"])
    # L_REC_0 is the list of recordings to be done BEFORE the simulations begin
    # L_REC is the list of recordings to be done DURING the simulations
    # L_REC_1 is the list of recordings to be done AFTER the simulations end
    #results is where all results will be stored

    global_var.current_timestep = 0
    # save the state of the system BEFORE the simulations begins
    for P_rec_i in L_REC_0:
        #P_rec_i contains all parameters for saving one single variable before the simulation
        single_save(Graphs, P_rec_i, results, global_var, internal_tick = 0)     # 0 means that the time step is before the simulations begin
        #writes data to              results

    for internal_tick in range(1,P_simulations["T"]+1):
        global_var.current_timestep = internal_tick

        if not global_var.stop_condition:
            # for each time step, advance the simulation one increment
            for P_rule in list_of_rules:
                single_update(Graphs, P_rule, global_var)
        
            # save the state of the system DURING the simulation
            # with the possibility to save the state every N steps via P_rec_i["DT"]
            for P_rec_i in L_REC:
                if(internal_tick%P_rec_i["DT"] == 0):
                    single_save(Graphs, P_rec_i, results, global_var, internal_tick = internal_tick)
                    #writes data to              results

        if global_var.stop_condition:
            for P_rec_i in L_REC:
                batch_save(Graphs, P_rec_i,results, global_var, internal_tick = internal_tick, T = P_simulations["T"])
                #writes lots of data to    results
            break

    # save the state of the system AFTER the simulations end
    for P_rec_i in L_REC_1:
        single_save(Graphs, P_rec_i, results, global_var, internal_tick = -1)     # -1 means that the time step is after the simulations end
        #writes data to              results


    convert_results_to_float(results)

    return results


def convert_results_to_float(results):
    """this function goes through all the attributes of results, and converts the arg.data into a np.float64 array"""
    attrs = [a for a in dir(results) if not a.startswith('_')]
    for attr in attrs:
        if(type(getattr(results, attr)) == np.ndarray):
            setattr(results, attr, getattr(results, attr).astype(getattr(results, attr)[0].dtype))



def single_update(G, P_rule, global_var = Global_Vars()):
    update_fct_name = P_rule["func"]
    try:
        update_fct_dict[update_fct_name](G, P_rule, global_var)  # Call the function
    except KeyError:
        print("I cannot update the following rule:" + update_fct_name)
        print("ERROR: single_update of rule(" + update_fct_name + ") NOT IMPLEMENTED YET, or something is broken in it")
    pass






