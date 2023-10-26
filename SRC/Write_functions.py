import numpy as np
import h5py

def add_dataset_to_group(simulation_group, name_dataset, values):
    if np.isscalar(values):
        simulation_group.create_dataset(name_dataset, data=values)
    # Check if the attribute value is a list of numbers (int or float)
    elif isinstance(values, list) and all(isinstance(i, (int, float)) for i in values):
        simulation_group.create_dataset(name_dataset, data=values)
    # Check if the attribute value is a list of lists of numbers (int or float)
    elif isinstance(values, list) and all(isinstance(i, list) for i in values):
        data_array = np.array(values)
        if len(data_array.shape) > 1:
            simulation_group.create_dataset(name_dataset, data=values)
    elif isinstance(values, np.ndarray):
        if len(values.shape) == 1:
            simulation_group.create_dataset(name_dataset, data=np.array(values, dtype=type(values[0])))
        elif len(values.shape) == 2:
            simulation_group.create_dataset(name_dataset, data=np.array(values, dtype=type(values[0][0])))
        else:
            # If the array has more than 2 dimensions, flatten it
            flattened_array = values.flatten()
            simulation_group.create_dataset(name_dataset, data=flattened_array.data)


def write_h5(i, P_rec, results):
    filename = P_rec["filename"]
    with h5py.File(filename, 'a') as hf:
        # Create a new simulation group within the file

        name_group = f"Simulation_{i}"
        if name_group in hf:
            del hf[name_group]
        simulation_group = hf.create_group(f"Simulation_{i}")
        # Iterate through the attributes of the object A
        for attr_name, attr_value in results.__dict__.items():
            data = np.array(attr_value.data.tolist())
            time = np.array(attr_value.time.tolist())
            #print("data = ", data, "\ntime = ", time)

            name_dataset_data = f"{attr_name}-data"
            name_dataset_time = f"{attr_name}-time"

            add_dataset_to_group(simulation_group, name_dataset_data, data)
            add_dataset_to_group(simulation_group, name_dataset_time, time)


# What follow will end up in a different file, with all the reading functions, and the one needed to parse and reconstruct the data. 

def load_data_from_h5(filename):
    data = {}
    with h5py.File(filename, 'r') as hf:
        for key in hf.keys():
            data[key] = {}
            for subkey in hf[key].keys():
                data[key][subkey] = np.array(hf[key][subkey])

    return data