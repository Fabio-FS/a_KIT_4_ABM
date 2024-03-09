import numpy as np
import igraph as ig



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
    #save_in_dict(P_rec, "DATA", RES)

def save_mean(G,P_rec):
    RES = np.mean(G[P_rec["layer"]].vs[P_rec["target"]])    
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def save_median(G,P_rec):
    RES = np.median(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def save_var(G,P_rec):
    RES = np.var(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

# add other polarization measures: Esteban Ray, std of pairwise differences, etc.

def save_max(G,P_rec):
    RES = np.max(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def save_min(G,P_rec):
    RES = np.min(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def save_frac(G,P_rec):
    RES =np.sum(np.array(G[P_rec["layer"]].vs[P_rec["target"]]) == P_rec["target_fraction"])/len(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def save_histogram(G,P_rec):
    RES = np.histogram(G[P_rec["layer"]].vs[P_rec["target"]])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def save_homophily(G,P_rec):
    RES = calc_H(G[P_rec["layer"]], P_rec["target"])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def save_fr_local(G,P_rec):
    RES = fr_local(G[P_rec["layer"]].vs[P_rec["target"]],P_rec["target"])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def save_graph(G,P_rec):
    # in case of several replicas, this overwrites the previous graph. I haven't looked into how to append a graph to a file yet.
    RES = ig.write_graph(G[0], P_rec["filename"])
    return RES
    #save_in_dict(P_rec, "DATA", RES)

def fr_local(g, name, value):
    # returns the fraction of neighbors with the same value  g is the graph  name is the name of the attribute  value is the value of the attribute
    # for each node, get the list of neighbors, and the fraction of neighbors with value = value, and store it in a list called RES

    RES = []
    
    for i in range(len(g.vs)):
        RES.append(np.array(g.vs[g.neighbors(i)][name]) == value)/len(g.neighbors(i))
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
    "fr_local" : save_fr_local,
    "graph" : save_graph
}



# writing functions for the hdf5 files
def save_scalar(name, i, n_trials, key, value):
    #print("this should not be used... only np arrays?")
    #print("key: " + key + " value: " + str(value))
    with h5py.File(name, 'a') as f:
        if key in f.keys():
            print("append")
            f[key][i] = value
        else:
            print("initialize")
            temp = np.empty(n_trials, dtype=type(value))
            temp[i] = value[0]
            print("temp=",temp)
            print("temp.shape=", temp.shape)
            print("temp.dtype=", temp.dtype)
            print("value=",value)
            print("i=",i)
            print("key = ", key)
        
            f.create_dataset(key, data=temp)
    return None

def write_array(name, key, i, value):
    print("saving_array")
    with h5py.File(name, 'a') as f:
        print("key = ", key)
        key = "Trial_" + str(i) + "_" + key
        f.create_dataset(key, data=value, chunks=True, maxshape=(None,))
        #f[key].resize((f[key].shape[0] + value.shape[0]), axis=0)
    #            f[key][-value.shape[0]:] = value
    #        else:
    #            
    #return None



#i, n_trials+1, P_rec, Data
def write_h5(i, n_trials, P_rec, results):

    print("parameters are: i=",i," n_trials=",n_trials," P_rec=",P_rec," results=",results)
    filename = P_rec["filename"]
    attributes =  [a for a in dir(results) if not a.startswith('_')]
    for name_attr in attributes:
        Obj = getattr(results,name_attr)
        L = len(Obj.data)
        print("L=",L)
        print(" I should enter in write array")

        if L > 1:
            # if the length of the array is > 1, save it as an array, and indipendently for each trial
            # the name of the array will be: T_$i$_name
            write_array(filename, name_attr, i, Obj)
            pass
        else:
            # if the length of the array is 1, create an array and save it inside the i-th trial elements of the array
            # the name of the array will be: name
            write_scalar(filename, name_attr, i, n_trials, Obj)
            pass


def write_scalar(filename, attribute_name, i, n_trials, Obj):
    """
    if "i" = 0, saves in the h5 file an array of length n_trials, where each entry is equal to Obj.data.
    if "i" > 0, saves in the h5 file the value Obj.data in the i-th entry of the array.
    """
    if i == 0:
        # if i = 0, initialize the array
        temp = np.empty(n_trials, dtype=type(Obj.data))
        temp[:] = Obj.data.copy()
        # save the array in the h5 file
        with h5py.File(filename, 'a') as f:
            f.create_dataset(attribute_name, data=temp)
    else:
        with h5py.File(filename, 'a') as f:
            f[attribute_name][i] = value

        


def save_scalar(name, i, n_trials, key, value):
    #print("this should not be used... only np arrays?")
    #print("key: " + key + " value: " + str(value))
    with h5py.File(name, 'a') as f:
        if key in f.keys():
            f[key][i] = value
        else:
            temp = np.empty(n_trials, dtype=type(value))
            temp[i] = value[0]
            print("temp=",temp)
            print("temp.shape=", temp.shape)
            print("temp.dtype=", temp.dtype)
            print("value=",value)
            print("i=",i)
            print("key = ", key)
        
            f.create_dataset(key, data=temp)
    return None







def tell_me_what_is_saved(name):
    with h5py.File(name, 'r') as f:
        for key in f.keys():
            print(key, f[key].shape)
    return None

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