def is_regular(g):
    #tests if a graph is regular
    #returns bool
    degrees = g.degree()
    return all(d == degrees[0] for d in degrees)

def row_mean_fast(arr):
    return np.mean(arr, axis = 1)

def row_mean_ragged(arr):
    return np.array([np.mean(row) for row in arr])

row_mean_fct = {1 : row_mean_fast, 0 : row_mean_ragged}