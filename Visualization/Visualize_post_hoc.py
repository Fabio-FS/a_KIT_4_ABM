import numpy as np
import igraph as ig
import matplotlib.pyplot as plt


def plot_scatterplot(graph, layer, variable_name, title, xlabel, ylabel, color = "black", size = 1, alpha = 0.4, save_as = "scatterplot.png"):
    # gets the list of the edges of the graph in the form of tuples
    edges = graph[layer].get_edgelist()
    # gets the list of the values of "variable_name" for each vertex in the tuples:

    V = graph[layer].vs[variable_name]
    #creates the list of the values of "variable_name" for each edge in the graph
    edge_values0 = [V[edge[0]] for edge in edges]
    edge_values1 = [V[edge[1]] for edge in edges]


    print("len edges0 = ", len(edge_values0))


    plt.scatter([edge_values0, edge_values1] , [edge_values1, edge_values0], c=color, s=size, alpha=alpha)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.savefig(save_as)
    plt.show()
    plt.close()

