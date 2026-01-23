import os
import importlib


local_path = os.path.dirname(os.path.abspath(__file__)) # Get the current script's directorys
folder_path = os.path.join(local_path, 'models') # Construct the folder_path by appending "/models" to the local_path
file_names = [f[:-3] for f in os.listdir(folder_path) if f.endswith('.py')] # Get a list of all files in the folder



update_fct_dict = {}
init_fct_dict = {}



# Import each module dynamically
for module_name in file_names:
    module = importlib.import_module(f'models.{module_name}')
    # If you want to import all names from the module, you can use vars() or globals()
    globals().update(vars(module))
    if hasattr(module, 'init_model') and callable(module.init_model):
        module.init_model(update_fct_dict, init_fct_dict)
        #init_model manipulates both dictionaries
        #and connects keys to actual functions, defined in subfolder models/