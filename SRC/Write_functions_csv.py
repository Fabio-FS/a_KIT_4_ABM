import numpy as np
import csv


def write_csv(P_rec, res):
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
        #print(fieldnames)
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        if not file_exists:
            writer.writeheader()

            row = {}
            for attr in fieldnames:
                row[attr] = getattr(getattr(res, attr), 'time')
            writer.writerow(row)

        row = {}
        for attr in fieldnames:
            row[attr] = getattr(getattr(res, attr), 'data')

        writer.writerow(row)