import numpy as np
import igraph as ig
import matplotlib.pyplot as plt

def load_F_of_t(namefile):
    filename = "namefile.txt"
    Xi = []
    XXi = []
    with open("fraction_of_infected.txt", "r") as file:
        for line in file:
            print(line)
            line = line.strip()  # Remove leading/trailing whitespace and newline characters
            if line.startswith("[") and line.endswith("]"):
                line = line[1:-1]  # Remove the brackets from the line
                elements = line.split(",\t")  # Split the line into elements based on ",\t"
                #print(elements)
                #print(type(elements))
                for element in elements:
                    if element.isdigit():  # Check if the element is a digit
                        Xi.append(float(element))
                    else:
                        Xi.append(float(element))
                        # Process non-numeric elements as needed
                        pass
            XXi.append(Xi)
    return np.array(XXi)