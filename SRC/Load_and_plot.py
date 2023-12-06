import numpy as np
import h5py


def load_data_from_h5(filename):
    data = {}
    with h5py.File(filename, 'r') as hf:
        for key in hf.keys():
            data[key] = {}
            for subkey in hf[key].keys():
                data[key][subkey] = np.array(hf[key][subkey])

    return data


def load_and_average(namefile, namevariable):
    A = load_data_from_h5(namefile)
    for i in range(len(A)):
        namef = "Simulation_" + str(i)
        named = namevariable + "-data"
        namet = namevariable + "-time" 
        if i == 0:
            data = np.array(A[namef][named])
            time = np.array(A[namef][namet])
        else:
            data += np.array(A[namef][named])
    return data/len(A), time


def load_end_states(namefile, namevariable):
    A = load_data_from_h5(namefile)
    for i in range(len(A)):
        namef = "Simulation_" + str(i)
        named = namevariable + "-data"
        if i == 0:
            data = np.array(A[namef][named][-1])
        else:
            data = np.vstack((data, np.array(A[namef][named][-1])))
    return data

def load_all(filepath, target):
    data = load_data_from_h5(filepath)
    RES = np.zeros(len(data))
    for i in range(len(data)):
        namef = "Simulation_" + str(i)
        named = target + "-data"
        RES[i] = np.array(data[namef][named][-1])
    return RES

def load_and_bin(filepath, target, bin_var, bins):
    V1 = load_all(filepath, target)
    V2 = load_all(filepath, bin_var)

    MEAN = np.zeros(len(bins)-1)
    MEDIAN = np.zeros(len(bins)-1)
    STD = np.zeros(len(bins)-1)
    
    for i in range(len(bins)-1):
        MEAN[i] = np.mean(V1[(V2>bins[i]) & (V2<bins[i+1])])
        MEDIAN[i] = np.median(V1[(V2>bins[i]) & (V2<bins[i+1])])
        STD[i] = np.std(V1[(V2>bins[i]) & (V2<bins[i+1])])
        
    return MEAN, MEDIAN, STD




### to be used for scalar only results

def load_res_scalars(filename):
    results = {}
    with h5py.File(filename, 'r') as file:
        for dataset_name in file:
            results[dataset_name] = file[dataset_name][:]

    return results



def load_all_scalars(filename, target):
    data = load_res_scalars(filename)
    target = "result_" + target
    return data[target]