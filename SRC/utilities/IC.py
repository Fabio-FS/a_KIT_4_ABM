import numpy as np
import igraph as ig
#from Clustering import *
from Hom_and_pol import *
import os


def set_disease_initial_condition(IC, attribute, g):
    if(IC["type"] == "random"):
        g.vs[attribute]=1
        # select N_pat_zero nodes at random and set them to 2:
        for gg in rn.sample(range(len(g.vs)), IC["N_pat_zero"]):

            g.vs[gg][attribute] = 2
    elif(IC["type"] == "central_geometric"):
        if IC["N_pat_zero"] != 1:
            print("\n")
            print("WARNING: Will always initialize with exactly 1 patient zero when type is central_geometric.")
            print("(Even though N_pat_zero is set to"+str(IC["N_pat_zero"])+")")
            print("\n")
        g.vs[attribute]=1
        # select N_pat_zero nodes at random and set them to 2:
        gg = int(np.floor(len(g.vs)/2))
        g.vs[gg][attribute] = 2
    elif(IC["type"] == "vector"):
        g.vs[attribute] = IC["values"]
    else:
        print("INITALIZATION RULE: " + IC["type"] + " NOT IMPLEMENTED YET! NUUUUUU")



# implementation of the initial conditions. Used for continuous dynamics

def set_continuous_initial_condition(IC, attribute, g, vector_from_init_fct = None):
    

    if not vector_from_init_fct is None:
        g.vs[attribute] = vector_from_init_fct
    elif(IC["distribution"] == "random_uniform"):
        #25-04-07: in older parameter files, this might still be "type"
        g.vs[attribute]=np.random.uniform(IC["Low"],IC["High"],len(g.vs))
    elif(IC["distribution"] == "random_beta"):
        # if alpha is not defined, sets it to be equal to beta, and viceversa
        if("alpha" not in IC):
            a, b = IC["beta"], IC["beta"]
        elif("beta" not in IC):
            a, b = IC["alpha"], IC["alpha"]
        else:
            a, b = IC["alpha"], IC["beta"]

        VECTOR_INIT = np.random.beta(a,b,len(g.vs))
        g.vs[attribute]= VECTOR_INIT
    elif(IC["dsitribution"] == "vector"):
        g.vs[attribute] = IC["values"]
    else:
        print("INITALIZATION RULE: " + IC["type"] + " NOT IMPLEMENTED YET! NUUUUUU")
        
    if(IC["homophily"]["Flag"] == True):
        # case 1 debug = False, return_H_hist = False
        dbg = IC["homophily"]["debug"]
        return_H_hist = IC["homophily"]["return_H_hist"]

        results = metropolis(   g, 
                                attribute =  attribute, 
                                hom_target = IC["homophily"]["hom_target"], 
                                N_steps = IC["homophily"]["steps"],
                                return_H_hist = return_H_hist,
                                dbg = dbg,
                                recenter = IC["homophily"]["recenter"])
        save_results_for_range_pol_hom(results, g, dbg, return_H_hist)



def save_results_for_range_pol_hom(results, g, debug, return_H_hist):
    if (not return_H_hist and not debug):
        rr = results
        pass
    if (return_H_hist and not debug):
        H_hist = results[1]
    if (not return_H_hist and debug):
        hist_B, A = results[1], results[2]
        pass
    if (return_H_hist and debug):
        H_hist = results[1]
        hist_B, A = results[2], results[3]

    if(return_H_hist):
        if os.path.exists("H_hist.csv"):
            os.remove("H_hist.csv")
        np.savetxt("H_hist.csv", H_hist, delimiter=",")
    
    if(debug):
        if os.path.exists("hist_B.csv"):
            os.remove("hist_B.csv")
        indx = calc_indxs(hist_B)
        hist_B2 = hist_B[indx,:]

        np.savetxt("hist_B.csv", hist_B2, delimiter=",")
        if os.path.exists("ADJ.csv"):
            os.remove("ADJ.csv")
        ig.Graph.write_adjacency(g, "ADJ.csv")


def calc_indxs(all_B):
    indx = []
    for i in range(int(np.log10(len(all_B)))):
        # append np.power(10,i) to indx
        indx.append(np.power(10,i))
        indx.append(2*np.power(10,i))
        indx.append(4*np.power(10,i))
    indx.append(len(all_B)-1)
    return indx