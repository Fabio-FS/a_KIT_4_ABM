import numpy as np
from Hom_and_pol import calc_homophily, calc_polarization

class Results():
    #an empty container for ALL results
    #saved values for each parameter will be an attribute of the Results object
    #e.g.  results.beta
    pass

class SingleVariable_Results:
    #each object of this class will store all results for a specific variable i
    #each object has two attributes: .time and .data
    def __init__(self,P_rec_i):
        empty_vector = np.empty(P_rec_i["total_count"], dtype=object)
        self.data = empty_vector.copy()
        self.time = P_rec_i["time_vector"]


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
    #where all results will be saved:
    results = Results()

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

        time_vector_i = []
        if(P_rec_i["END"] == "True"):
            L_REC_1.append(P_rec_i)
            P_rec_i["total_count"] += 1
            time_vector_i.append(-1)
        if(P_rec_i["BEGIN"] == "True"):
            L_REC_0.append(P_rec_i)
            P_rec_i["total_count"] += 1
            time_vector_i.append(0)
        if(P_rec_i["DT"] > 0):
            L_REC.append(P_rec_i)
            P_rec_i["total_count"] += int(T_max/P_rec_i["DT"])
            time_vector_i.extend(   np.arange(1,int(T_max/P_rec_i["DT"])+1,P_rec_i["DT"]).tolist())
        
        P_rec_i["time_vector"] = np.array(time_vector_i)
        
        if P_rec_i["total_count"] > 0:
            # if at least one of the conditions is satisfied, initialize the results for single variable
            setattr(results, P_rec_i["column_name"], SingleVariable_Results(P_rec_i))

    return L_REC_0, L_REC, L_REC_1, results      # list of recordings at the beginning of the simulation, during the simulation, and at the end of the simulation



def single_save(G, P_rec_i, results, internal_tick = -10):
    saving_fct_name = P_rec_i["func"]
    try:
        saving_function = saving_dictionary[saving_fct_name]    
    except KeyError:
        print("ERROR: " + saving_fct_name + " not found")
    
    RES = saving_function(G, P_rec_i)
    col_name  = P_rec_i["column_name"]

    idx = (internal_tick  +   (P_rec_i["END"]=="True"))   *  (1 - (internal_tick == -1))
    getattr(results,col_name).data[idx] = RES             # RES  is saved in position idx of results.name.data

def batch_save(G, P_rec_i, results, internal_tick = -10, T=500):
    saving_fct_name = P_rec_i["func"]
    try:
        saving_function = saving_dictionary[saving_fct_name]    
    except KeyError:
        print("ERROR: " + saving_fct_name + " not found")
    
    RES = saving_function(G, P_rec_i)
    col_name  = P_rec_i["column_name"]
    
    next_internal_tick = (internal_tick // P_rec_i["DT"])*P_rec_i["DT"] + P_rec_i["DT"]
    idx = (next_internal_tick  +   (P_rec_i["END"]=="True"))
    
    if type(RES) == list:
        getattr(results,col_name).data[idx:] = [RES for _ in np.arange(next_internal_tick,T+1,P_rec_i["DT"])]            # RES  is saved in all subsequent positions of results.name.data
        if P_rec_i["END"] == "True":
            getattr(results,col_name).data[0] = [RES]
    else:    
        getattr(results,col_name).data[idx:] = RES             # RES  is saved in all subsequent positions of results.name.data
        if P_rec_i["END"] == "True":
            getattr(results,col_name).data[0] = RES


def save_ALL(G,P_rec):
    RES = G[P_rec["layer"]].vs[P_rec["attribute"]]
    return RES

def save_mean(G,P_rec):
    RES = float(np.mean(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

def save_median(G,P_rec):
    RES = float(np.median(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

def save_var(G,P_rec):
    RES = float(np.var(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

# add other polarization measures: Esteban Ray, std of pairwise differences, etc.

def save_max(G,P_rec):
    RES = float(np.max(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

def save_min(G,P_rec):
    RES = float(np.min(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

def save_frac(G,P_rec):
    RES = float(np.sum(np.array(G[P_rec["layer"]].vs[P_rec["attribute"]]) == P_rec["attribute_value"])/len(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

def save_histogram(G,P_rec):
    RES = np.histogram(G[P_rec["layer"]].vs[P_rec["attribute"]])
    return RES

def save_homophily(G,P_rec):
    RES = float(calc_homophily(G[P_rec["layer"]], P_rec["attribute"]))
    return RES

def save_homophily_rescaled(G,P_rec):
    RES = float(calc_homophily(G[P_rec["layer"]], P_rec["attribute"], flag = 1))
    return RES

def save_homophily_non_rescaled(G,P_rec):
    RES = float(calc_homophily(G[P_rec["layer"]], P_rec["attribute"], flag = 2))
    return RES


def save_fr_local(G,P_rec):
    RES = fr_local(G[P_rec["layer"]].vs[P_rec["attribute"]],P_rec["attribute"])
    return RES

def fr_local(g, name, value):
    # returns the fraction of neighbors with the same value. g is the graph,  name is the name of the attribute,  value is the value of the attribute
    RES = []
    
    for i in range(len(g.vs)):
        RES.append(np.array(g.vs[g.neighbors(i)][name]) == value)/len(g.neighbors(i))
    return RES

def save_pol(G,P_rec):
    RES = calc_polarization(G[P_rec["layer"]], P_rec["attribute"])
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
    "homophily_rescaled" : save_homophily_rescaled,
    "homophily_non_rescaled" : save_homophily_non_rescaled,
    "hom" : save_homophily,
    "fr_local" : save_fr_local,
    "pol" : save_pol,
    "polarization" : save_pol
}
