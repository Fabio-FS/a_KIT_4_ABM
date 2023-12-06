global_path = 'C:/Users/nel_t/Documents/WORK/PROJECTS/a_KIT_4_ABM'   # <--- CHANGE THIS PATH to the one on your computer
code_path   = global_path + '/SRC'

#global_path = 'C:/Users/nel_t/Documents/WORK/PROJECTS/a_KIT_4_ABM'   # <--- CHANGE THIS PATH to the one on your computer
#code_path   = global_path + '/SRC'


# import the KIT_4_ABM package
import sys
sys.path.append(code_path)
import KIT_4_ABM as kit
import numpy as np
import argparse

 
def reset_param(P_lay_original,P_sim_original,P_dyn_original,P_rec_original,ab):
    P_lay = P_lay_original.copy()
    P_dyn = P_dyn_original.copy()
    P_sim = P_sim_original.copy()
    P_rec = P_rec_original.copy()
    P_dyn["Dynamic_0"]["BEHAVIOR"]["IC"]["alpha"] = ab
    return P_lay, P_sim, P_dyn, P_rec


def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='My Simulation')

    # Add command-line arguments for A and B
    parser.add_argument('--va', type=float, help='Vaccines availability')
    parser.add_argument('--s2vh', type=float, help='max P of vaccination')
    parser.add_argument('--wsp', type=float, help='rewiring probability of WS')
    
    # Parse the command-line arguments
    args = parser.parse_args()

    va_ff = args.va             # Vaccines availability
    S2V_h_ff = args.s2vh       # max P of vaccination
    WSP_ff = args.wsp           # WS probability of rewiring

    np.random.seed(42)

    P_lay, P_dyn, P_sim, P_rec = kit.import_parameters("PAR_SIRV.json")
    P_lay["Layer_0"]["P"] = WSP_ff


    N = 3
    d = 0.25/N
    
    range_pol = np.linspace(0+d,0.25-d, N-1)
    value_ab = 1/(8*range_pol)-0.5

    P_rec["filename"] = P_rec["filename"]+ "_VA=" + str(va_ff)+ "_S2Vh="
    P_rec["filename"] = P_rec["filename"] + str(S2V_h_ff) + "_WSP=" + str(WSP_ff) + ".h5"

    P_dyn["Dynamic_0"]["VA"] = va_ff

    P_lay_original = P_lay.copy()
    P_dyn_original = P_dyn.copy()
    P_sim_original = P_sim.copy()
    P_rec_original = P_rec.copy()

    n_trials = 1
    count = 0
    for i in range(n_trials):
        for ab in value_ab:
            print("Trial: " + str(count) + " of " + str(n_trials*N), "alpha = " + str(ab))
            
            P_lay, P_sim, P_dyn, P_rec = reset_param(P_lay_original, P_sim_original, P_dyn_original, P_rec_original, ab)

            res = kit.run_sim(P_lay, P_dyn, P_sim, P_rec)
            kit.write_h5(count, P_rec,res)
            #kit.write_h5(count, P_rec, res)
            #if P_sim["save_settings"] == "Round":
            #    pass
            #else:
            #it.save_res_scalars(count, res, P_rec["filename"], n_trials*len(value_ab)) 
            count += 1
        #if P_sim["save_settings"] == "Round":
            # rename the res
            #print(type(res))
            #kit.add_to_res(res, res2, str(round(ab,4)))
            #kit.write_h5(count, P_rec, res2) 
            #count += 1



if __name__ == '__main__':
    main()
