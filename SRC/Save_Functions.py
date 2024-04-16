import numpy as np
import igraph as ig
#from Clustering import calc_H
from Hom_and_pol import calc_homophily, calc_polarization


class Saves:
    pass



#  ██████  ███████  ██████  ██████  ██████  ██████  ██ ███    ██  ██████  ███████ 
#  ██   ██ ██      ██      ██    ██ ██   ██ ██   ██ ██ ████   ██ ██       ██      
#  ██████  █████   ██      ██    ██ ██████  ██   ██ ██ ██ ██  ██ ██   ███ ███████ 
#  ██   ██ ██      ██      ██    ██ ██   ██ ██   ██ ██ ██  ██ ██ ██    ██      ██ 
#  ██   ██ ███████  ██████  ██████  ██   ██ ██████  ██ ██   ████  ██████  ███████ 


def init_recordings(P_recordings, T_max):

    # initialize the recordings
    L_REC_0 = []                        # list of recordings at the beginning of the simulation
    L_REC   = []                        # list of recordings during the simulation
    L_REC_1 = []                        # list of recordings at the end of the simulation

    for i in range(P_recordings["N"]):
        # for each recording i, read the parameters and initialize the recording
        P_rec_i = P_recordings["Recording_" + str(i)]
        # for each of recording, I calculate how often it appears:
        P_rec_i["total_count"] = 0

        if "BEGIN" not in P_rec_i:
            # if BEGIN is not specified, set it to False
            P_rec_i["BEGIN"] = "False"
        if "END" not in P_rec_i:
            # if END is not specified, set it to False
            P_rec_i["END"] = "False"
        if "DT" not in P_rec_i:
            # if DT is not specified, set it to 0
            P_rec_i["DT"] = 0
            # I add to count the number of times 1:DT:T_max is divisible by DT

        
        #print("P_rec_i[totla_count]=", P_rec_i["total_count"])
        if(P_rec_i["BEGIN"] == "True"):
            L_REC_0.append(P_rec_i)
            P_rec_i["total_count"] += 1
        #print("P_rec_i[totla_count]=", P_rec_i["total_count"])
        if(P_rec_i["END"] == "True"):
            L_REC_1.append(P_rec_i)
            P_rec_i["total_count"] += 1
        #print("P_rec_i[totla_count]=", P_rec_i["total_count"])
        if(P_rec_i["DT"] > 0):
            L_REC.append(P_rec_i)
            P_rec_i["total_count"] += int(T_max/P_rec_i["DT"])
        #print("P_rec_i[total_count]=", P_rec_i["total_count"])
        
        if P_rec_i["total_count"] > 0:
            # if at least one of the conditions is satisfied, set the count to 0
            P_rec_i["count"] = 0

    return L_REC_0, L_REC, L_REC_1      # list of recordings at the beginning of the simulation, during the simulation, and at the end of the simulation



#single_save(G, record, internal_tick, Delta_T = record["DT"])
def single_save(G, P_rec, results, internal_tick = -10):
    name = P_rec["func"]
    #print("count = ", P_rec["count"], "before single save is executed")
    try:
        function = saving_dictionary[name]    
    except KeyError:
        print("ERROR: " + name + " NOT IMPLEMENTED YET! NUUUUUU")
    
    RES = function(G, P_rec)
    name  = P_rec["name"]
        
    if P_rec["count"] == 0:
        # In Data, I create the field P_rec["name"], with length L = P_rec["total_count"]. Each element of Data.P_rec["name"] is a copy of RES
        empty_vector = np.empty(P_rec["total_count"], dtype=object)
        dummy = Saves()
        dummy.data = empty_vector.copy()
        dummy.time = empty_vector.copy()

        setattr(results, name, dummy)
        #setattr(results, name, empty_vector.copy())


    
    getattr(results,name).data[P_rec["count"]] = RES             # RES  is saved in position P_rec["count"] of results.name.data,
    getattr(results,name).time[P_rec["count"]] = internal_tick   # TIME is saved in position P_rec["count"] of results.name.time

    #print("name_var=",  name_var,  "var =", RES)
    #print("name_time=", name_time, "tic =", internal_tick)

    #print(getattr(Data,name_var))
    #print(getattr(Data,name_time))
    
    P_rec["count"] += 1                                         # P_rec["count"] is increased to save in the right position in the next iteration



# measure for the health dynamic: 
#   - hight of the peak of the wave
#   - time above a threshold
#   - 




def save_ALL(G,P_rec):
    RES = G[P_rec["layer"]].vs[P_rec["target"]]
    return RES

def save_mean(G,P_rec):
    RES = np.mean(G[P_rec["layer"]].vs[P_rec["target"]])    
    return RES

def save_median(G,P_rec):
    RES = np.median(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES

def save_var(G,P_rec):
    RES = np.var(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES

# add other polarization measures: Esteban Ray, std of pairwise differences, etc.

def save_max(G,P_rec):
    RES = np.max(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES

def save_min(G,P_rec):
    RES = np.min(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES

def save_frac(G,P_rec):
    RES =np.sum(np.array(G[P_rec["layer"]].vs[P_rec["target"]]) == P_rec["target_fraction"])/len(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES

def save_histogram(G,P_rec):
    RES = np.histogram(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES

def save_homophily(G,P_rec):
    RES = calc_homophily(G[P_rec["layer"]], P_rec["target"])
    return RES

def save_fr_local(G,P_rec):
    RES = fr_local(G[P_rec["layer"]].vs[P_rec["target"]],P_rec["target"])
    return RES

def fr_local(g, name, value):
    # returns the fraction of neighbors with the same value. g is the graph,  name is the name of the attribute,  value is the value of the attribute
    RES = []
    
    for i in range(len(g.vs)):
        RES.append(np.array(g.vs[g.neighbors(i)][name]) == value)/len(g.neighbors(i))
    return RES

def save_pol(G,P_rec):
    RES = calc_polarization(G[P_rec["layer"]], P_rec["target"])
    return RES

saving_dictionary = {
    "ALL" : save_ALL,
    "avg" : save_mean,
    "mean" : save_mean,
    "median" : save_median,
    "var" : save_var,
    "variance" : save_var,
    "max" : save_max,
    "maximum" : save_max,
    "min" : save_min,
    "minimum" : save_min,
    "frac" : save_frac,
    "fraction" : save_frac,
    "hist" : save_histogram,
    "histogram" : save_histogram,
    "homophily" : save_homophily,
    "hom" : save_homophily,
    "fr_local" : save_fr_local,
    "pol" : save_pol,
    "polarization" : save_pol
}