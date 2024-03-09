import numpy as np
import os
import importlib


local_path = os.path.dirname(os.path.abspath(__file__)) # Get the current script's directorys
folder_path = os.path.join(local_path, 'models') # Construct the folder_path by appending "/models" to the local_path
file_names = [f[:-3] for f in os.listdir(folder_path) if f.endswith('.py')] # Get a list of all files in the folder




update_dictionary = {}
init_single_rule = {}



# Import each module dynamically
for module_name in file_names:
    module = importlib.import_module(f'models.{module_name}')
    # If you want to import all names from the module, you can use vars() or globals()
    globals().update(vars(module))
    if hasattr(module, 'init_model') and callable(module.init_model):
        module.init_model(update_dictionary, init_single_rule)




# Ask Arne how to make the code expandable without expanding the code. I want users to be able to add new rules without having to modify the code.
# Update schedule?
    # how to do the savings in case of evolutionary game theory? (several waves of infection)

# Same models with different communication regimes (121, 12M, M21)




# Ability to read csv files of graphs
# Temporal network 
# Ability to add new save functions



