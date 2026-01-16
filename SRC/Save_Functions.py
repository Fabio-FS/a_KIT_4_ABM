import numpy as np
from Hom_and_pol import calc_homophily, calc_polarization
from Write_functions import *

#   ████    ██         ██      ████     ████    ██████    ████   
#  ██  ██   ██        ████    ██  ██   ██  ██   ██       ██  ██  
#  ██       ██       ██  ██   ██       ██       ██       ██      
#  ██       ██       ██████    ████     ████    ████      ████   
#  ██       ██       ██  ██       ██       ██   ██           ██  
#  ██  ██   ██       ██  ██   ██  ██   ██  ██   ██       ██  ██  
#   ████    ██████   ██  ██    ████     ████    ██████    ████   


class Results():
    #an empty container for ALL results
    #saved values for each parameter will be an attribute of the Results object
    #e.g.  results.beta
    pass

class SingleVariable_Results:
    #each object of this class will store all results for a specific variable i
    #each object has two attributes: .time and .data
    def __init__(self,P_rec_i):
        empty_vector = np.empty(P_rec_i["total_count"], dtype=np.float64)
        self.data = empty_vector.copy()
        self.time = P_rec_i["time_vector"]


#  ██████  ███████  ██████  ██████  ██████  ██████  ██ ███    ██  ██████  ███████ 
#  ██   ██ ██      ██      ██    ██ ██   ██ ██   ██ ██ ████   ██ ██       ██      
#  ██████  █████   ██      ██    ██ ██████  ██   ██ ██ ██ ██  ██ ██   ███ ███████ 
#  ██   ██ ██      ██      ██    ██ ██   ██ ██   ██ ██ ██  ██ ██ ██    ██      ██ 
#  ██   ██ ███████  ██████  ██████  ██   ██ ██████  ██ ██   ████  ██████  ███████ 


def init_recording(P_record, T_max):

    # initialize the recordings
    L_REC_0 = []                        # what to record this at the BEGINNING of the simulation
    L_REC   = []                        # what to record this DURING the simulation
    L_REC_1 = []                        # what to record at the END of the simulation
    #where all results will be saved:
    results = Results()

    for i in range(P_record["N"]):
        # for each recording i, read the parameters and initialize the recording
        P_rec_i = P_record["Recording_" + str(i)]
        # for each of recording, I calculate how often it appears:
        P_rec_i["total_count"] = 0

        if "BEGIN" not in P_rec_i:
            # if BEGIN is not specified, set it to False
            P_rec_i["BEGIN"] = False
        if "END" not in P_rec_i:
            # if END is not specified, set it to False
            P_rec_i["END"] = False
        if "DT" not in P_rec_i:
            # if DT is not specified, set it to 0
            P_rec_i["DT"] = 0
            # I add to count the number of times 1:DT:T_max is divisible by DT

        time_vector_i = []
        if(P_rec_i["END"] == True):
            L_REC_1.append(P_rec_i)
            P_rec_i["total_count"] += 1
            time_vector_i.append(-1)
        if(P_rec_i["BEGIN"] == True):
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



def single_save(G, P_rec_i, results, global_var, internal_tick = -10):
    #writes values to the Results object
    #the Results object is later exported to csv
    calc_fct_name = P_rec_i["func"]
    try:
        statistic_function = saving_dictionary[calc_fct_name]    
    except KeyError:
        print("ERROR: " + calc_fct_name + " not found")
    
    RES = statistic_function(G, P_rec_i, global_var)
    col_name  = P_rec_i["column_name"]

    idx = (internal_tick  +   (P_rec_i["END"]==True))   *  (1 - (internal_tick == -1))
    getattr(results,col_name).data[idx] = RES             # RES  is saved in position idx of results.name.data

def batch_save(G, P_rec_i, results, global_var, internal_tick = -10, T=500):
    #writes many values at once to the Results object
    #the Results object is later exported to csv
    calc_fct_name = P_rec_i["func"]
    try:
        statistic_function = saving_dictionary[calc_fct_name]    
    except KeyError:
        print("ERROR: " + calc_fct_name + " not found")
    
    RES = statistic_function(G, P_rec_i, global_var)
    col_name  = P_rec_i["column_name"]
    
    next_internal_tick = (internal_tick // P_rec_i["DT"])*P_rec_i["DT"] + P_rec_i["DT"]
    idx = (next_internal_tick  +   (P_rec_i["END"]==True))
    
    if type(RES) == list:
        getattr(results,col_name).data[idx:] = [RES for _ in np.arange(next_internal_tick,T+1,P_rec_i["DT"])]            # RES  is saved in all subsequent positions of results.name.data
        if P_rec_i["END"] == True:
            getattr(results,col_name).data[0] = [RES]
    else:    
        getattr(results,col_name).data[idx:] = RES             # RES  is saved in all subsequent positions of results.name.data
        if P_rec_i["END"] == True:
            getattr(results,col_name).data[0] = RES

                                                                                                                                                                                    
#   ████    ██████     ██     ██████    ████     ████    ██████    ████     ████     ████  
#  ██  ██     ██      ████      ██       ██     ██  ██     ██       ██     ██  ██   ██  ██ 
#  ██         ██     ██  ██     ██       ██     ██         ██       ██     ██       ██     
#   ████      ██     ██████     ██       ██      ████      ██       ██     ██        ████  
#      ██     ██     ██  ██     ██       ██         ██     ██       ██     ██           ██ 
#  ██  ██     ██     ██  ██     ██       ██     ██  ██     ██       ██     ██  ██   ██  ██ 
#   ████      ██     ██  ██     ██      ████     ████      ██      ████     ████     ████  


#  ██████   ██  ██   ██  ██    ████    ██████    ████     ████    ██  ██    ████   
#  ██       ██  ██   ███ ██   ██  ██     ██       ██     ██  ██   ███ ██   ██  ██  
#  ██       ██  ██   ██████   ██         ██       ██     ██  ██   ██████   ██      
#  ████     ██  ██   ██████   ██         ██       ██     ██  ██   ██████    ████   
#  ██       ██  ██   ██ ███   ██         ██       ██     ██  ██   ██ ███       ██  
#  ██       ██  ██   ██  ██   ██  ██     ██       ██     ██  ██   ██  ██   ██  ██  
#  ██        ████    ██  ██    ████      ██      ████     ████    ██  ██    ████  



def pass_ALL(G,P_rec, global_var = None):
    if P_rec.get("source",None) == "global_var":
        return getattr(global_var,P_rec["attribute"])
    RES = G[P_rec["layer"]].vs[P_rec["attribute"]]
    return RES

def pass_ALL_as_int(G,P_rec, global_var = None):
    RES = np.array(G[P_rec["layer"]].vs[P_rec["attribute"]]).astype(int).tolist()
    return RES

def calc_mean(G,P_rec, global_var = None):
    if P_rec.get("source",None) == "global_var":
        return float(np.mean(getattr(global_var,P_rec["attribute"])))
    return float(np.mean(G[P_rec["layer"]].vs[P_rec["attribute"]]))

def calc_mean_subgroup(G,P_rec,global_var = None):
    subgroup_idcs = np.where(    np.array(G[P_rec["layer"]].vs[P_rec["attribute_subgroup_1"]])    ==     P_rec["attribute_subgroup_1_value"])[0]
    if P_rec.get("source",None) == "global_var":
        with np.errstate(invalid="ignore"):
            return float(np.mean(getattr(global_var,P_rec["attribute"])[subgroup_idcs]))
    with np.errstate(invalid="ignore"):
        return float(np.mean(G[P_rec["layer"]].vs[subgroup_idcs][P_rec["attribute"]] ))

    

def calc_median(G,P_rec, global_var = None):
    if P_rec.get("source",None) == "global_var":
        return float(np.median(getattr(global_var,P_rec["attribute"])))
    RES = float(np.median(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

def calc_var(G,P_rec, global_var = None):
    RES = float(np.var(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

# add other polarization measures: Esteban Ray, std of pairwise differences, etc.

def calc_max(G,P_rec, global_var = None):
    RES = float(np.max(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

def calc_min(G,P_rec, global_var = None):
    RES = float(np.min(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    return RES

def calc_frac(G,P_rec, global_var = None):
    if P_rec.get("source",None) == "global_var":
        return float(   np.mean(getattr(global_var,P_rec["attribute"])    ==    P_rec["attribute_value"])    )
    return float(   np.mean(np.array(G[P_rec["layer"]].vs[P_rec["attribute"]])   ==   P_rec["attribute_value"])   )
    #RES = float(np.sum(np.array(G[P_rec["layer"]].vs[P_rec["attribute"]]) == P_rec["attribute_value"])/len(G[P_rec["layer"]].vs[P_rec["attribute"]]))
    #return RES

def calc_frac_subgroup(G, P_rec, global_var = None):

    not_all_conditions_applied = True
    subgroup_index = np.array(list(range(G[P_rec["layer"]].vcount())))
    all_idcs = np.array(list(range(G[P_rec["layer"]].vcount())))
    true_counter = 0
    condition_idx = 1   
    while not_all_conditions_applied:
        true_counter += 1
        try:
            new_condition_idcs = all_idcs[np.array(G[P_rec["layer"]].vs[P_rec[f"attribute_subgroup_{condition_idx}"]]) == P_rec[f"attribute_subgroup_{condition_idx}_value"]]
            subgroup_index = np.intersect1d(subgroup_index, new_condition_idcs )
        except:
            not_all_conditions_applied = False
        condition_idx += 1

    if len(subgroup_index) == 0:
        return 0
    
    subgroup_size = len(subgroup_index)/G[P_rec["layer"]].vcount()
    
    return float( np.mean(np.array(G[P_rec["layer"]].vs[subgroup_index][P_rec["attribute"]])   ==   P_rec["attribute_value"]) * subgroup_size  )

def calc_homophily_wrapper(G,P_rec, global_var = None):
    RES = float(calc_homophily(G[P_rec["layer"]], P_rec["attribute"], recenter_flag = P_rec["recenter_flag"], qs = P_rec["qs"], is_category = P_rec["is_category"]))
    return RES

def calc_homophily_recentered_wrapper(G,P_rec, global_var = None):
    RES = float(calc_homophily(G[P_rec["layer"]], P_rec["attribute"], recenter_flag = 1, qs = P_rec["qs"], is_category = P_rec["is_category"]))
    return RES

def calc_homophily_non_recentered_wrapper(G,P_rec, global_var = None):
    RES = float(calc_homophily(G[P_rec["layer"]], P_rec["attribute"], recenter_flag = 2, qs = P_rec["qs"], is_category = P_rec["is_category"]))
    return RES


def calc_fr_local(G,P_rec, global_var = None):
    RES = fr_local(G[P_rec["layer"]].vs[P_rec["attribute"]],P_rec["attribute"])
    return RES

def fr_local(g, attribute, value, global_var = None):
    # returns the fraction of neighbors with the same value. g is the graph,  name is the name of the attribute,  value is the value of the attribute
    RES = []
    
    for i in range(len(g.vs)):
        RES.append(np.array(g.vs[g.neighbors(i)][attribute]) == value)/len(g.neighbors(i))
    return RES

def calc_pol(G,P_rec, global_var = None):
    RES = calc_polarization(G[P_rec["layer"]], P_rec["attribute"])
    return RES

saving_dictionary = {
    "ALL" : pass_ALL,
    "ALL_int" : pass_ALL_as_int,
    "avg" : calc_mean,
    "mean" : calc_mean,
    "mean_subgroup" : calc_mean_subgroup,
    "median" : calc_median,
    "var" : calc_var,
    "variance" : calc_var,
    "max" : calc_max,
    "maximum" : calc_max,
    "min" : calc_min,
    "minimum" : calc_min,
    "frac" : calc_frac,
    "fraction" : calc_frac,
    "frac_subgroup" : calc_frac_subgroup,
    "hom" : calc_homophily_wrapper, 
    "homophily" : calc_homophily_wrapper,
    "homophily_recentered" : calc_homophily_recentered_wrapper,
    "homophily_non_recentered" : calc_homophily_non_recentered_wrapper,
    "fr_local" : calc_fr_local,
    "pol" : calc_pol,
    "polarization" : calc_pol
}
