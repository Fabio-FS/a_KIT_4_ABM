#!/usr/bin/env python
# coding: utf-8


import numpy as np
import igraph as ig
from Init_and_Dynamics import *
from Save_Functions import *
from Write_functions_csv import *
import sys

import json

class Global_Vars:
    def __init__(self):
        self.stop_condition = False


def import_parameters(namefile):
    # namefile is a json file in the JSON format.
    # first I import everything into data, then I split data into P_lay, P_dyn, P_sim, P_rec

    with open(namefile) as f:
        data = json.load(f)

    P_lay = data["P_layer"]
    P_dyn = data["P_dynamic"]
    P_sim = data["P_simulations"]
    P_rec = data["P_recordings"]

    return P_lay, P_dyn, P_sim, P_rec

def run_sim(P_layer, P_dynamic, P_simulations, P_recordings, return_G = False):

    # initialize the graph creating all the needed layers. for each layer create the network
    G = init_graph(P_layer) # G is a list of graphs.

    global_var = Global_Vars()
    #global_var contains variables that all functions should have access to in principle
    #for the moment it's only stop_condition

    # initializes the dynamic on the graph and returns a list of rules for the updating function.
    list_of_rules = init_rules(G, P_dynamic,global_var)
    
    Data = run_temporal_evolution(G, list_of_rules, P_simulations, P_recordings,global_var)

    if return_G:
        return Data, G
    return Data

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
            g = ig.Graph.Lattice(dim=[p_lay_i["Lx"], p_lay_i["Ly"]], circular=False, nei = p_lay_i["nei"])
        elif(p_lay_i["type"] == "Moore_Lattice"):
            #a lattice with Moore neighborhood, i.e. 8 neighbors
            #number of agents
            N = p_lay_i["Lx"]*p_lay_i["Ly"]
            #set up adjacency matrix
            A = np.zeros(shape = (N,N))

            #horizontal: +1
            #vertical:   +x  (conflict with horiz. if x=1)
            #NW-SE:      +x+1
            #NE-SW:      +x-1 (conflict with horiz. if x=2)

            if p_lay_i["Lx"] == 2:
                #special treatment because of NE-SW connections
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
                    else:
                        A[i,2*(i//2)-2:2*(i//2)-2+6] = np.array([1,1,1,0,1,1])
                if p_lay_i["circular"] == "True":
                    A[N-1,1] = 1    #connnecting last agent with 2nd
                    A[N-2,0] = 1    #connecting 2nd to last agent with first
                    A[N-1,0] = 1    #connnecting last agent with 1st
                    A[N-2,1] = 1    #connecting 2nd to last agent with 2nd
                A += A.T
                A = np.clip(A,a_min=0,a_max=1)
            else:   #xdim is not 2
                A += np.diag(np.ones(shape = N-1),1 )    #horizontal conections, one off the diagonal
                A[-p_lay_i["Lx"]+np.arange(p_lay_i["Lx"], N, 1)   ,   np.arange(p_lay_i["Lx"], N, 1)]   =   1     #vertical
                A[-p_lay_i["Lx"]-1 + np.arange(p_lay_i["Lx"]+1,N,1)   ,   np.arange(p_lay_i["Lx"]+1,N,1)]  =  1    #NW-SE
                A[-p_lay_i["Lx"]+1 + np.arange(p_lay_i["Lx"],N,1)   ,   np.arange(p_lay_i["Lx"],N,1)]  =  1    #NE-SW
                
                A += A.T

                if p_lay_i["Lx"]>1:   #otherwise I will delete the vertical connections here, because for x=2:   x-1 = 1
                    A[np.arange(p_lay_i["Lx"]-1,N-1,p_lay_i["Lx"]),np.arange(p_lay_i["Lx"],N,p_lay_i["Lx"])] = 0 # horizontal correction
                A[-p_lay_i["Lx"]-1+np.arange(2*p_lay_i["Lx"],N,p_lay_i["Lx"] )   ,   np.arange(2*p_lay_i["Lx"],N,p_lay_i["Lx"] )] = 0   #NW-SE correction
                A[np.arange(0,N,p_lay_i["Lx"])   ,   p_lay_i["Lx"]-1+np.arange(0,N,p_lay_i["Lx"])] = 0 #NW-SE   correction

                A *= A.T

                if p_lay_i["circular"] == "True":
                    #A_pre = np.copy(A[:,:])
                    #print(A)
                    A[np.arange(0,p_lay_i["Lx"],1)     ,    np.arange(N-p_lay_i["Lx"],N,1)]    =   1    #adding vertical connections
                    #print(np.arange(0,p_lay_i["Lx"],1)     ,    np.arange(N-p_lay_i["Lx"],N,1))
                    #print(A)
                    #print(A-A_pre)
                    if p_lay_i["Lx"] > 1:
                        A[np.arange(p_lay_i["Lx"]-1,N,p_lay_i["Lx"])   ,   np.arange(0,N,p_lay_i["Lx"])]   =   1   #adding horizontal connections
                        #if x dimension is bigger than 1, otherwise I introduce self loops here
                        A[np.arange(p_lay_i["Lx"]-1,N-p_lay_i["Lx"],p_lay_i["Lx"])    ,   1+np.arange(p_lay_i["Lx"]-1,N-p_lay_i["Lx"],p_lay_i["Lx"])]   =   1   #adding NW-SE connections on the right edge
                        A[N-1,0] = 1   #connecting right-bottom to top-left
                        A[p_lay_i["Lx"]-1,N-p_lay_i["Lx"]] = 1
                        #print(p_lay_i["Lx"]-1,N-p_lay_i["Lx"])
                        #print(np.arange(2*p_lay_i["Lx"]-1,N,p_lay_i["Lx"])    ,   np.arange(0,N-p_lay_i["Lx"]+1,p_lay_i["Lx"]))
                        A[np.arange(2*p_lay_i["Lx"]-1,N,p_lay_i["Lx"])    ,   np.arange(0,N-p_lay_i["Lx"],p_lay_i["Lx"])]   =   1   #adding NE-SW connections on the right edge
                        A[np.arange(N-p_lay_i["Lx"],N-1,1)    ,    np.arange(1,p_lay_i["Lx"],1)]   =  1   #adding NE-SW connecs in bottom row
                        A[np.arange(N-p_lay_i["Lx"]+1,N,1)    ,    np.arange(0,p_lay_i["Lx"]-1,1)]   =  1   #adding NW-SE connecs in bottom row
                        #print(np.arange(N-p_lay_i["Lx"],N-1,1)    ,    np.arange(1,p_lay_i["Lx"],1))
                        #print(np.arange(N-p_lay_i["Lx"]+1,N,1)    ,    np.arange(0,p_lay_i["Lx"]-1,1))
                    A += A.T
                    A = np.clip(A,a_min=0,a_max=1)
                    #print(A)
            A = np.asarray(A, dtype=int)
            g = ig.Graph.Adjacency(A, mode="undirected")
        elif(p_lay_i["type"] == "WS"):
            if(np.power(p_lay_i["L"],p_lay_i["D"]) != N):
                print("GRAPH SIZE WARNING: L^D != N: " + str(np.power(p_lay_i["L"],p_lay_i["D"])) + " != " + str(N))
            g = ig.Graph.Watts_Strogatz(dim=p_lay_i["D"], size=p_lay_i["L"], nei=p_lay_i["NFN"], p=p_lay_i["P"])
        elif p_lay_i["type"] == "WS_MS":
            #Watts-Strogatz network with Maslov-Sneppen-like rewiring
            if(np.power(p_lay_i["L"],p_lay_i["D"]) != N):
                print("GRAPH SIZE WARNING: L^D != N: " + str(np.power(p_lay_i["L"],p_lay_i["D"])) + " != " + str(N))
            g = ig.Graph.Watts_Strogatz(dim=p_lay_i["D"], size=p_lay_i["L"], nei=p_lay_i["NFN"], p=0)
        elif(p_lay_i["type"] == "kRRG"):
            g = ig.Graph.K_Regular(n = N, k = p_lay_i["k"])
        elif(p_lay_i["type"] == "2islands"):

            p = p_lay_i["p"]
            k = p_lay_i["k"]
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
            print("GRAPH: " + p_lay_i["type"] + " not implemented yet")

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

def run_temporal_evolution(G, list_of_rules, P_simulations, P_recordings,global_var):
    # G is the list of graph-layers, each item is one graph
    # list_of_rules is a list of dictionaries, each containing parameters for one updating function
    # P_recordings is the dictionary with the parameters for the recordings

    np.set_printoptions(threshold=sys.maxsize)
    #necessary for networks with >999 agents.
    #otherwise numpy prints every array as [y_0, y_1, ..., y_n]

    L_REC_0, L_REC, L_REC_1, results = init_recordings(P_recordings, P_simulations["T"])
    # L_REC_0 is the list of recordings to be done BEFORE the simulations begin
    # L_REC is the list of recordings to be done DURING the simulations
    # L_REC_1 is the list of recordings to be done AFTER the simulations end
    #results is where all results will be stored

    # save the state of the system BEFORE the simulations begins
    for P_rec_i in L_REC_0:
        #P_rec_i contains all parameters for saving one single variable before the simulation
        single_save(G, P_rec_i, results, internal_tick = 0)     # 0 means that the time step is before the simulations begin

    for internal_tick in range(1,P_simulations["T"]+1):

        if not global_var.stop_condition:
            # for each time step, advance the simulation one increment
            for P_rule in list_of_rules:
                single_update(G, P_rule, global_var)
        
            # save the state of the system DURING the simulation
            # with the possibility to save the state every N steps via P_rec_i["DT"]
            for P_rec_i in L_REC:
                if(internal_tick%P_rec_i["DT"] == 0):
                    single_save(G, P_rec_i, results, internal_tick = internal_tick)

        if global_var.stop_condition:
            for P_rec_i in L_REC:
                batch_save(G, P_rec_i,results, internal_tick = internal_tick, T = P_simulations["T"])
            break

    # save the state of the system AFTER the simulations end
    for P_rec_i in L_REC_1:
        single_save(G, P_rec_i, results, internal_tick = -1)     # -1 means that the time step is after the simulations end


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






