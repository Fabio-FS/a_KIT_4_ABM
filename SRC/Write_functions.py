import csv
import json
import numpy as np
import h5py
import KIT_4_ABM as kit

def write_csv(P_rec, res, sweep_parameters = None):

    filename = P_rec["filename"] + ".csv"
    try:
        with open(filename, 'r') as existing_file:
            reader = csv.reader(existing_file)
            file_exists = any(reader)
    except FileNotFoundError:
        file_exists = False
        
    with open(filename, 'a', newline='') as csvfile:
        # Extract attribute names ending with '.data'
        fieldnames = [attr for attr in dir(res) if hasattr(getattr(res, attr), 'data')]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if not file_exists:
            writer.writeheader()

            row = {}
            for attr in fieldnames:
                row[attr] = getattr(getattr(res, attr), 'time')
            writer.writerow(row)

            if not sweep_parameters is None:
                for attr in fieldnames:
                    row[attr] = sweep_parameters
                writer.writerow(row)


        row = {}
        for attr in fieldnames:
            row[attr] = getattr(getattr(res, attr), 'data')
            
        writer.writerow(row)


def sanity_check_csv(P_rec, res, sweep_parameters = None, trim = 20):

    filename = P_rec["filename"] + ".csv"
    try:
        with open(filename, 'r') as existing_file:
            reader = csv.reader(existing_file)
            file_exists = any(reader)
    except FileNotFoundError:
        file_exists = False
        
    with open(filename, 'a', newline='') as csvfile:
        # Extract attribute names ending with '.data'
        fieldnames = [attr for attr in dir(res) if hasattr(getattr(res, attr), 'data')]

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if not file_exists:
            writer.writeheader()

            row = {}
            for attr in fieldnames:
                attr_time = getattr(getattr(res, attr), 'time')
                if len(attr_time) > trim:
                    attr_time = attr_time[:trim]
                row[attr] = attr_time
            writer.writerow(row)

            if not sweep_parameters is None:
                for attr in fieldnames:
                    row[attr] = sweep_parameters
                writer.writerow(row)


        row = {}
        for attr in fieldnames:
            attr_data = getattr(getattr(res, attr), 'data')
            if not sweep_parameters is None:
                if len(attr_data.shape) > 1:
                    if attr_data.shape[1] > trim:
                        attr_data = attr_data[:,:trim]
                #if instead, the shape has only 1 dim, these are just one value per parameter combination
                #and it shouldn't be trimmed
            elif attr_data.shape[0] > trim:
                attr_data = attr_data[:trim]
            row[attr] = attr_data
            
        writer.writerow(row)

def create_h5(P_record):

    filename = P_record["filename"] + ".h5"

    with h5py.File(filename, "w", libver="latest") as f:

        if P_record.get("live_monitoring",True):
            f.swmr_mode = True
        #needed for live monitoring.
        #in swmr mode, multiple processes can read 

        #create a dataset for sweep_codes
        total_sweep_steps = 1
        if not P_record.get("sweep_codes") is None:
            f.create_dataset(
                "sweep_codes",
                data = P_record.get("sweep_codes")
            )
            total_sweep_steps = P_record.get("sweep_codes").shape[1]
            P_record["total_sweep_steps"] = total_sweep_steps

        for i in range(P_record["N"]):
            P_rec_i = P_record["Recording_" + str(i)]

            grp = f.require_group(P_rec_i["column_name"])

            time_vector = P_rec_i["time_vector"]

            grp.create_dataset(
                "timestamps",
                data = time_vector
            )

            N_values = P_rec_i["N_values"]

            grp.create_dataset(
                "values",
                shape=   (total_sweep_steps, 0   , len(time_vector), N_values),      # initial shape : (parameter configs, trials, timesteps, values saved per config-trial-timestep)
                maxshape=(total_sweep_steps, None, len(time_vector), N_values),      # unlimited 2nd dim
                chunks=  (total_sweep_steps, 1   , len(time_vector), N_values),      # one trial per chunk
                dtype="float64"
            )

def append_h5(P_rec, res):

    filename = P_rec["filename"] + ".h5"

    with h5py.File(filename, "a") as f:
        fieldnames = [attr for attr in dir(res) if hasattr(getattr(res, attr), 'data')]

        for attr in fieldnames:
            values = f[f"{attr}/values"]
            values.resize(values.shape[1] + 1, axis=1)
            values[:,-1,:,:] = getattr(getattr(res, attr), 'data').reshape(values.shape[0],values.shape[2],values.shape[3])
            #[parameter combination, trial, timestamps, N_values per network]
            f.flush() #prevents corruption in case of aborting simulation with Ctrl C
        