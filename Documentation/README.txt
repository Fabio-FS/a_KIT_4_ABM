 ▄▄▄▄▄▄▄▄▄▄▄       ▄    ▄  ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄▄       ▄         ▄       ▄▄▄▄▄▄▄▄▄▄▄  ▄▄▄▄▄▄▄▄▄▄   ▄▄       ▄▄ 
▐░░░░░░░░░░░▌     ▐░▌  ▐░▌▐░░░░░░░░░░░▌▐░░░░░░░░░░░▌     ▐░▌       ▐░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░▌ ▐░░▌     ▐░░▌
▐░█▀▀▀▀▀▀▀█░▌     ▐░▌ ▐░▌  ▀▀▀▀█░█▀▀▀▀  ▀▀▀▀█░█▀▀▀▀      ▐░▌       ▐░▌     ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░▌░▌   ▐░▐░▌
▐░▌       ▐░▌     ▐░▌▐░▌       ▐░▌          ▐░▌          ▐░▌       ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌▐░▌ ▐░▌▐░▌
▐░█▄▄▄▄▄▄▄█░▌     ▐░▌░▌        ▐░▌          ▐░▌          ▐░█▄▄▄▄▄▄▄█░▌     ▐░█▄▄▄▄▄▄▄█░▌▐░█▄▄▄▄▄▄▄█░▌▐░▌ ▐░▐░▌ ▐░▌
▐░░░░░░░░░░░▌     ▐░░▌         ▐░▌          ▐░▌          ▐░░░░░░░░░░░▌     ▐░░░░░░░░░░░▌▐░░░░░░░░░░▌ ▐░▌  ▐░▌  ▐░▌
▐░█▀▀▀▀▀▀▀█░▌     ▐░▌░▌        ▐░▌          ▐░▌           ▀▀▀▀▀▀▀▀▀█░▌     ▐░█▀▀▀▀▀▀▀█░▌▐░█▀▀▀▀▀▀▀█░▌▐░▌   ▀   ▐░▌
▐░▌       ▐░▌     ▐░▌▐░▌       ▐░▌          ▐░▌                    ▐░▌     ▐░▌       ▐░▌▐░▌       ▐░▌▐░▌       ▐░▌
▐░▌       ▐░▌     ▐░▌ ▐░▌  ▄▄▄▄█░█▄▄▄▄      ▐░▌                    ▐░▌     ▐░▌       ▐░▌▐░█▄▄▄▄▄▄▄█░▌▐░▌       ▐░▌
▐░▌       ▐░▌     ▐░▌  ▐░▌▐░░░░░░░░░░░▌     ▐░▌                    ▐░▌     ▐░▌       ▐░▌▐░░░░░░░░░░▌ ▐░▌       ▐░▌
 ▀         ▀       ▀    ▀  ▀▀▀▀▀▀▀▀▀▀▀       ▀                      ▀       ▀         ▀  ▀▀▀▀▀▀▀▀▀▀   ▀         ▀ 
                                                                                                                  


-------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------------------



 ██████  ██    ██ ███████ ██████  ██    ██ ██ ███████ ██     ██ 
██    ██ ██    ██ ██      ██   ██ ██    ██ ██ ██      ██     ██ 
██    ██ ██    ██ █████   ██████  ██    ██ ██ █████   ██  █  ██ 
██    ██  ██  ██  ██      ██   ██  ██  ██  ██ ██      ██ ███ ██ 
 ██████    ████   ███████ ██   ██   ████   ██ ███████  ███ ███  
                                                                
                                                                


This envirorment allows you to simulate several dynamics on a pseudo-multiplex network. Each dynamic can be defined
on a specific layer of the multiplex network iteself.

All the parameters for the simulation must be defined in one parameter file.
STEP 1: write a parameter.json file (see PAR_list) for a comprehensive list of implemented dynamic, savings and networks
STEP 2: load the parameter file using
    P_layers, P_dynamics, P_simulation, P_recordings = import_parameters("parameter.json")
STEP 2.5: adjust the parameters as needed.
STEP 3: run the simulations using
    run_sim(P_layers, P_dynamics, P_simulation, P_recordings)
STEP 4: read the results in the txt files generated while runnign the simulations

------------------------------------------- 
-------- STEP 1 --------------------------- 
------------------------------------------- 

The parameter file is a json file. It MUST contains 4 dictionaries:
    1) P_simulations    - contains the details of the simulation
                            - number of timesteps
                            - number of waves if "evo-SIR" (not implemented yet)
    2) P_layer          - contains info about the number of layers in the multiplex network, and all the details to generate them
    3) P_dynamic        - contains info about the number of dynamics, and the layer where they are living in. And the name of the attribute they are reading from other dynamics
    4) P_recordings     - contains info about what to save. specifically "func" tells the name of the function to be called, and "namefile" the name of the file where everything will be saved on

------------------------------------------- 
-------- STEP 2 --------------------------- 
------------------------------------------- 

Calling the function import_parameters will return 4 dictionaries with the needed parameters to run the simulations.
This is not embedded directly within the run_sim, to allow the user to change dynamically one of the parameters, without modifying the parameter file.

------------------------------------------- 
-------- STEP 3 --------------------------- 
------------------------------------------- 

Calling run_sim will do 3 things:

1) Build a pseudo-multiplex network (uses info from P_layer)
2) Imprint the dynamics in the PMN (uses info from P_dynamic)
3) Run the dynamics (for a number of step and waves specified in P_simulation) and save on the fly the outcome of the simulation (using info from P_recordings)

-----
# G = init_graph(P_layer) # G is a list of graphs.

Since in igraph there is no native way to filter edges by attribute having several edges corresponding to different layers of the multiplex network, is rather slow.
For this reason I implemented each layer of the network as independent graph, and saved them in the list of graphs G.


-----
# list_of_rules, list_of_layers = imprint_rules(G, P_dynamic)

Since we are now working with several stacked networks, each dynamic must specify on which layer it's operating. "home_layer".
In case of dynamics that read parameters defined in a different dynamic, as for the SIR-beta, the parameter file must specify:
    - in which layer to find the "behavior" variable.
    - the name of the behavior variable.
It's important to double check that the names matches. Since there is no internal consistency check, and I have no clue what the outcome would be otherwise.
In principle it's possible to have several dynamics on the same layers. as long as the names of the variables the dynamic is operating with are different that is not a problem. (shouldn't be a problem)
    print("The graph is now completed. now running the simulations...")

-----
# run_temporal_evolution(G, list_of_rules, list_of_layers, P_simulations, P_recordings)
selfexplanatory for now.

