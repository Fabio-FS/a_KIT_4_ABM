import numpy as np
import igraph as ig
from Clustering import calc_H
import h5py







#  ██████  ███████  ██████  ██████  ██████  ██████  ██ ███    ██  ██████  ███████ 
#  ██   ██ ██      ██      ██    ██ ██   ██ ██   ██ ██ ████   ██ ██       ██      
#  ██████  █████   ██      ██    ██ ██████  ██   ██ ██ ██ ██  ██ ██   ███ ███████ 
#  ██   ██ ██      ██      ██    ██ ██   ██ ██   ██ ██ ██  ██ ██ ██    ██      ██ 
#  ██   ██ ███████  ██████  ██████  ██   ██ ██████  ██ ██   ████  ██████  ███████ 


def init_recordings(P_recordings):

    # initialize the recordings
    L_REC_0 = []                        # list of recordings at the beginning of the simulation
    L_REC   = []                        # list of recordings during the simulation
    L_REC_1 = []                        # list of recordings at the end of the simulation

    for i in range(P_recordings["N"]):
        # for each recording i, read the parameters and initialize the recording
        P_rec_i = P_recordings["Recording_" + str(i)]
        if(P_rec_i["BEGIN"] == "True"):
            L_REC_0.append(P_rec_i)
        if(P_rec_i["END"] == "True"):
            L_REC_1.append(P_rec_i)
        if(P_rec_i["DT"] > 0):
            L_REC.append(P_rec_i)

    return L_REC_0, L_REC, L_REC_1      # list of recordings at the beginning of the simulation, during the simulation, and at the end of the simulation



#single_save(G, record, internal_tick, Delta_T = record["DT"])
def single_save(G, P_rec, tick, Delta_T = -10):
    name = P_rec["func"]
    try:
        function = saving_dictionary[name]
        function(G, P_rec)  # Call the function and save the results in P_rec["DATA"]
    except KeyError:
        print("ERROR: " + name + " NOT IMPLEMENTED YET! NUUUUUU")
    pass




def save_in_dict(temp, key, value):
    if key in temp:
    # if it exists, check if "X" is a list
        if isinstance(temp[key], list):
            # if it is a list, append the values
            temp[key].append(value)
        # if it is a vector, append the values
        elif isinstance(temp[key], np.ndarray):
            temp[key] = np.append(temp[key], value)
        # if it is a scalar, change it to a list and append the values
        else:
            temp[key] = [temp[key], value]
    # if the key "X" does not exist, create it and append the values
    else:
        temp[key] = value





def save_ALL(G,P_rec):
    RES = G[P_rec["layer"]].vs[P_rec["target"]]
    save_in_dict(P_rec, "DATA", RES)

def save_mean(G,P_rec):
    RES = np.mean(G[P_rec["layer"]].vs[P_rec["target"]])    
    save_in_dict(P_rec, "DATA", RES)

def save_median(G,P_rec):
    RES = np.median(G[P_rec["layer"]].vs[P_rec["target"]])
    save_in_dict(P_rec, "DATA", RES)

def save_var(G,P_rec):
    RES = np.var(G[P_rec["layer"]].vs[P_rec["target"]])
    save_in_dict(P_rec, "DATA", RES)

def save_max(G,P_rec):
    RES = np.max(G[P_rec["layer"]].vs[P_rec["target"]])
    save_in_dict(P_rec, "DATA", RES)

def save_min(G,P_rec):
    RES = np.min(G[P_rec["layer"]].vs[P_rec["target"]])
    save_in_dict(P_rec, "DATA", RES)

def save_frac(G,P_rec):
    RES =np.sum(np.array(G[P_rec["layer"]].vs[P_rec["target"]]) == P_rec["target_fraction"])/len(G[P_rec["layer"]].vs[P_rec["target"]])
    save_in_dict(P_rec, "DATA", RES)

def save_histogram(G,P_rec):
    RES = np.histogram(G[P_rec["layer"]].vs[P_rec["target"]])
    save_in_dict(P_rec, "DATA", RES)

def save_homophily(G,P_rec):
    RES = calc_H(G[P_rec["layer"]], P_rec["target"])
    save_in_dict(P_rec, "DATA", RES)

def save_fr_local(G,P_rec):
    RES = fr_local(G[P_rec["layer"]].vs[P_rec["target"]])
    save_in_dict(P_rec, "DATA", RES)

def save_graph(G,P_rec):
    # in case of several replicas, this overwrites the previous graph. I haven't looked into how to append a graph to a file yet.
    RES = ig.write_graph(G[0], P_rec["filename"])
    save_in_dict(P_rec, "DATA", RES)

def fr_local(g, name, value):
    # returns the fraction of neighbors with the same value
    # g is the graph
    # name is the name of the attribute
    # value is the value of the attribute
    # for each node, get the list of neighbors, and the fraction of neighbors with value = value, and store it in a list called RES

    RES = []
    
    for i in range(len(g.vs)):
        RES.append(np.array(g.vs[g.neighbors(i)][name]) == value)/len(g.neighbors(i))
    return RES




saving_dictionary = {
    "ALL" : save_ALL,
    "mean" : save_mean,
    "median" : save_median,
    "var" : save_var,
    "max" : save_max,
    "min" : save_min,
    "frac" : save_frac,
    "histogram" : save_histogram,
    "homophily" : save_homophily,
    "fr_local" : save_fr_local,
    "graph" : save_graph
}









def save_scalar(name, i, n_trials, key, value):
    with h5py.File(name, 'a') as f:
        if np.isscalar(value):
            if key in f.keys():
                f[key][i] = value
            else:
                temp = np.zeros(n_trials)
                temp[i] = value
                f.create_dataset(key, data=temp)
    return None

def save_array(name, i, key, value):
    with h5py.File(name, 'a') as f:
        key = "Trial_" + str(i) + "_" + key
        if isinstance(value, np.ndarray):
            if key in f.keys():
                f[key].resize((f[key].shape[0] + value.shape[0]), axis=0)
                f[key][-value.shape[0]:] = value
            else:
                f.create_dataset(key, data=value, chunks=True, maxshape=(None,))
    return None


def save_all(i, n_trials, name, **kwargs):
    for key, value in kwargs.items():
        if isinstance(value, np.ndarray):
            save_array(name, i, key, value)
        elif np.isscalar(value):
            save_scalar(name, i, n_trials, key, value)



def load_array(name, key):
    """
    loads a dataset from a hdf5 file
    """
    with h5py.File(name, 'r') as f:
        if key in f.keys():
            return f[key][:]
        else:
            print("this dataset does not exist")
            return None