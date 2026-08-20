import numpy as np
import igraph as ig


def n_simple_rewire(n_edges,p, size = 1):
    # calculates what to put for   n   in g.rewire(n = ?, mode = "simple")   as to be consistent with the Watts-Strogatz algorithm
    # the function g.rewire affects exactly and determinstically    2*n    edges.

    # However, in the Watts Strogatz model,  for   N    edges    and a rewiring probability   p   , on average p*N edges will be affected.
    # The exact value in each instance comes from the probability distribution    binomial(N,p)
    # The following formula captures the stochastic nature and is close to half the average of the above distribution.
    # close means: off by less than 0.25  (so the total number of affected edges is off by less than 0.5)

    #N:     number of edges
    #p:     rewiring probability
    #size:  how many values to generate

    return (    0.5* np.random.binomial ( n = n_edges, p=p, size = size )            ).astype(int)

def MooreLattice_adjacency(Lx,Ly,circular):
    #a lattice with Moore neighborhood, i.e. 8 neighbors
    #a Cellular automaton
    #parameters are "Lx", "Ly", "circular"

    #number of agents
    N = Lx*Ly
    #set up adjacency matrix
    A = np.zeros(shape = (N,N))
    
    #horizontal: +1
    #vertical:   +x  (conflict with horiz. if x=1)
    #NW-SE:      +x+1
    #NE-SW:      +x-1 (conflict with horiz. if x=2)
    
    #most examples given below are for the case N=9, Lx=3

    #→→→→ horizontal conections +1
    A[ np.where( np.arange(0,N) % Lx != Lx-1) [0]   ,   np.where( np.arange(0,N) % Lx != Lx-1) [0] + 1 ]   =   1
    #e.g.      0,1,3,4,6,7                                                     1,2,4,5,7,8
    #   0,1,2,3,4,5,6,7,8   --(%Lx != Lx-1)-->    0,1,3,4,6,7   --(add +1)-->  1,2,4,5,7,8

    #↓↓↓↓  vertical    +x
    A[np.arange(0, N-Lx,   1)   ,   np.arange(Lx,   N, 1)]   =   1
    #e.g.    0,1,2,3,4,5                  3,4,5,6,7,8
    
    # ↘↘↘↘  NW-SE       +x+1
    A[ np.where( np.arange(0,N-Lx) % Lx != Lx-1) [0]   ,   np.where( np.arange(0,N-Lx) % Lx != Lx-1)[0] + Lx+1 ]   =   1
    #e.g.    0,1,3,4                    4,5,7,8
    
    # ↙↙↙↙ NE-SW       +x-1
    A[ np.where( np.arange(0,N-Lx) % Lx != 0) [0]   ,   np.where( np.arange(0,N-Lx) % Lx != 0) [0] + Lx-1 ]   =   1
    #e.g.    1,2,4,5                  3,4,6,7
    
    #making adjacency matrix symmetric
    A += A.T
    
    #adding connections for circular lattices
    if circular == True:
        
        #1. ↓↓↓↓ adding vertical connections between first and last line
        A[np.arange(0,Lx,1)     ,    np.arange(N-Lx,N,1)]    =   1
        #e.g.       0,1,2                      6,7,8
        #
        #2. ↘↘↘↘ adding NW-SE connections on the right edge
        A[np.arange(Lx-1,N,Lx)    ,   np.mod(np.arange(Lx,N+1,Lx),N)]   =   1
        #e.g.       [2,5,8]                            [3,6,(9mod9=0)]
        #
        #3. ↙↙↙↙ adding NE-SW connections on right edge
        A[np.arange(Lx-1,N,Lx)    ,   np.mod(np.arange(-Lx,N-Lx,Lx),N)]   =   1
        #e.g.        2,5,8                           6,0,3   (-3,0,3)
        #
        #4. ↘↘↘↘ adding NW-SE connections in bottom row, except bottom-right
        A[np.arange(N-Lx,N-1,1)    ,    np.arange(1,Lx,1)]   =  1
        #e.g.        6,7                      1,2
        #5. ↙↙↙↙ adding NW-SE connecs in bottom row, except bottom-left
        A[np.arange(N-Lx+1,N,1)    ,    np.arange(0,Lx-1,1)]   =  1
        #e.g.        7,8                           0,1
        
        #6.→→→→ horizontal conections on right edge   -Lx+1;   only if Lx>1
        A[np.where( np.arange(0,N) % Lx == Lx-1)[0], np.where( np.arange(0,N) % Lx == Lx-1) [0] -Lx+1] = Lx > 1
        #if Lx == 1   this would add self-neighboring nodes:    0,1,2,3,4,5,6,7     --(%Lx == Lx-1)-->   0,1,2,3,4,5,6,7   --(add -Lx+1)-->   0,1,2,3,4,5,6,7
        #if Lx == 2   these nodes are already connected:        0,1,2,3,4,5,6,7     --(%Lx == Lx-1)-->   1,3,5,7           --(add -Lx+1)-->   0,2,4,6
        #if Lx >= 3   connect edges horizontally                0,1,2,3,4,5,6,7,8   --(%Lx == Lx-1)-->   2,5,8             --(add -Lx+1)-->   0,3,6
        #in other words: Lx==1 ->  right edge is left edge
        #                Lx==2 ->  right edge and left edge already connected
        #                Lx==3 ->  right edge and left edge not yet connected

        A += A.T
        A = np.clip(A,a_min=0,a_max=1)
    A = np.asarray(A, dtype=int)
    return A


def MS_rewire(g, g0, dim, p_rewire, size = (10,10), max_steps = 3000):
    #this function rewires a graph randomly.
    #while preserving the degree distribution
    #this is achieved by using a Maslov-Sneppen-inspired rewiring algorithm
    #Borrowing from Watts-Strogatz, we define a fixed and a loose end for each edge.
    #When rewiring, the fixed end of edge A will be wired to the loose end of edge B and vice versa
    #Fixed and loose ends are for now only determined for regular graphs: 1) circle networks and 2) Lattices with von Neumann or Moore neighborhood

    N = g.vcount()
    k = np.mean(g.degree())

    #effective probability of rewiring           (in WS, rewiring to the same node is possible)
    p_effective = p_rewire   *  ( (N-k-2) / (N-k-1) )

    #how many edges should actually be rewired?
    #WS is a stochastic process. Each edge will be rewired with probability p_effective
    #the total number of edges rewired follows a binomial distribution
    n_edges_rewired = (    np.random.binomial ( n = g.ecount(), p= p_effective, size = 1 )            ).astype(int)[0]
    n_edges_remaining = g.ecount() - n_edges_rewired

    #determining which end of each edge is fixed and which is open to rewiring
    edges_unordered = np.array([e.tuple for e in g.es])
    edges = np.zeros(edges_unordered.shape,dtype = int)

    for i,e in enumerate(edges_unordered):
        if dim == 1:
            e_fix   =  min(e)  if  (min(e) - (max(e) - N))   >   (max(e) - min(e))  else  max(e)
            e_loose =  min(e)  if  e_fix != min(e)                                  else  max(e)
        elif dim == 2:
            if (e[0]//size[0]) == (e[1]//size[0] + 1)%size[1]:
                # e[0] is one line below e[1]
                e_fix = e[1]
                e_loose = e[0]
            elif (e[1]//size[0]) == (e[0]//size[0] + 1)%size[1]:
                # e[1] is one line below e[0]
                e_fix = e[0]
                e_loose = e[1]
            elif e[0]%size[1] == (e[1] + 1)%size[1]:
                #e[0] is to the right of e[1]
                e_fix = e[1]
                e_loose = e[0]                        
            else:
                e_fix = e[0]
                e_loose = e[1]
        else:
            print("Rewiring not yet designed for graphs with dimensions above 2.")

        edges[i] = np.array([e_fix,e_loose])

    n_steps = 0
    n_overlap = g.ecount()
    while n_overlap > n_edges_remaining and n_steps < max_steps:
        idcs = np.random.choice(g.ecount(), size = 2, replace = False)

        if edges[idcs[0]][0] != edges[idcs[1]][1] and edges[idcs[1]][0] != edges[idcs[0]][1]:
            #no self loops

            nghbrs0 = g.neighbors(edges[idcs[0]][0])
            nghbrs1 = g.neighbors(edges[idcs[1]][0])

            if edges[idcs[1]][1] not in nghbrs0     or     ((edges[idcs[1]][1] in nghbrs0) and edges[idcs[1]][1] == edges[idcs[0]][1]) :
                #new link is not neighbor                    new link IS neighbor          but       new link is equal to old link
                #new node cannot be an existing neighbor, unless we cut the ties to it in the same step
                if edges[idcs[0]][1] not in nghbrs1     or     ((edges[idcs[0]][1] in nghbrs1) and edges[idcs[1]][1] == edges[idcs[0]][1]) :
                    #new link is not neighbor                    new link IS neighbor          but       new link is equal to old link
                    #new node cannot be an existing neighbor, unless we cut the ties to it in the same step

                    g.delete_edges([tuple(edges[idcs[0]]),tuple(edges[idcs[1]])])

                    g.add_edges([(edges[idcs[0]][0],edges[idcs[1]][1]),(edges[idcs[1]][0],edges[idcs[0]][1] )])
                    n_overlap = ig.intersection([g,g0]).ecount()

                    edges[idcs[0]][1], edges[idcs[1]][1] = edges[idcs[1]][1], edges[idcs[0]][1]
        n_steps += 1
    return g



# ██  ██   ██████   ██████   ██   ██   ████    █████    ██  ██  
# ███ ██   ██         ██     ██   ██  ██  ██   ██  ██   ██ ██   
# ██████   ██         ██     ██   ██  ██  ██   ██  ██   ████    
# ██████   ████       ██     ██ █ ██  ██  ██   █████    ███     
# ██ ███   ██         ██     ███████  ██  ██   ████     ████    
# ██  ██   ██         ██     ███ ███  ██  ██   ██ ██    ██ ██   
# ██  ██   ██████     ██     ██   ██   ████    ██  ██   ██  ██  


def init_graph(P_net):
    G = []
    for i in range(P_net["N_layers"]):
        P_layer = P_net[f"Layer_{i}"]
        N = P_net["N_nodes"]

        if(P_layer["type"] == "ER_m"):
            g = ig.Graph.Erdos_Renyi(n=N, p=P_layer["link_m"]/(N-1))
        elif(P_layer["type"] == "ER_p"):
            g = ig.Graph.Erdos_Renyi(n=N, p=P_layer["p"])
        elif(P_layer["type"] == "file"):
            #if P_layer has not the field "format" try to load the file with igraph.Load without forcing the format:
            if("format" not in P_layer):
                g = ig.Graph.Load(P_layer["filename"])
            else:
                g = ig.Graph.Read_Ncol(P_layer["filename"], format = P_layer["format"])
        elif(P_layer["type"] == "Lattice"):
            if(P_layer["Lx"]*P_layer["Ly"] != N):
                print("GRAPH SIZE WARNING: Lx*Ly != N: " + str(P_layer["Lx"]*P_layer["Ly"]) + " != " + str(N)) 
            g = ig.Graph.Lattice(dim=[P_layer["Lx"], P_layer["Ly"]], circular=P_layer.get("circular",False), nei = P_layer.get("nei",1))
        elif(P_layer["type"] == "Moore_Lattice"):
            #a lattice with Moore neighborhood, i.e. 8 neighbors
            #a Cellular automaton
            #parameters are "Lx", "Ly", "circular"

            A = MooreLattice_adjacency(P_layer["Lx"],P_layer["Ly"],P_layer.get("circular",False))
            g = ig.Graph.Adjacency(A, mode="undirected")

        elif(P_layer["type"] == "Moore_Lattice_WS_MS"):
            #a lattice with Moore neighborhood, i.e. 8 neighbors, then rewired Maslov-Sneppen-like
            #a Cellular automaton
            #parameters are "Lx", "Ly", "circular"

            A = MooreLattice_adjacency(P_layer["Lx"],P_layer["Ly"],P_layer.get("circular",False))
            g = ig.Graph.Adjacency(A, mode="undirected")
            g0 = ig.Graph.Adjacency(A, mode="undirected")

            g = MS_rewire(g,g0,dim = 2, p_rewire = P_layer.get("p_rewire",0), size = (P_layer["Lx"],P_layer["Ly"]))

        elif(P_layer["type"] == "WS"):
            if(np.power(P_layer["L"],P_layer["D"]) != N):
                print("GRAPH SIZE WARNING: L^D != N: " + str(np.power(P_layer["L"],P_layer["D"])) + " != " + str(N))
            g = ig.Graph.Watts_Strogatz(dim=P_layer["D"], size=P_layer["L"], nei=P_layer["n_neighbors"], p=P_layer["p_rewire"])
        elif(P_layer["type"] == "kRRG"):
            g = ig.Graph.K_Regular(n = N, k = P_layer["k"])
        elif(P_layer["type"] == "BA"):
            g = ig.Graph.Barabasi(n = N, m =  P_layer["m"])
        elif(P_layer["type"] == "2islands"):

            p = P_layer["p"]
            k = P_layer["k"]
            #island1 = ig.Graph.Erdos_Renyi(n=N//2, p=p_lay_i["k"]/(N//2-1))
            #island2 = ig.Graph.Erdos_Renyi(n=N-N//2, p=p_lay_i["k"]/((N-N//2)-1))

            membership = np.empty(N)
            membership[:N//2] = np.full(N//2,fill_value=0)
            membership[N//2:] = np.full(N-N//2,fill_value=1)

            gs0 = N//2
            gs1 = N - N//2
            success = False
            while success == False:

                #print(2*p*k/(N-1), 2*(1-p)*k/(N-1), p*k/(N-1) + (1-p)*k/(N-1))
                #print(np.log(N)/N)
                adj_matrix00 = np.random.choice([0,1],size = (gs0,gs0), p=[1-2*(1-p)*k/(N-1), 2*(1-p)*k/(N-1)])
                adj_matrix11 = np.random.choice([0,1],size = (gs1,gs1), p=[1-2*(1-p)*k/(N-1), 2*(1-p)*k/(N-1)])
                adj_matrix01 = np.random.choice([0,1],size = (gs0,gs1), p=[1-2*p*k/(N-1), 2*p*k/(N-1)])

                adj_matrix00 = np.tril(adj_matrix00, -1) + np.tril(adj_matrix00, -1).T
                adj_matrix11 = np.tril(adj_matrix11, -1) + np.tril(adj_matrix11, -1).T
                #adj_matrix01 = np.tril(adj_matrix01, -1) + np.tril(adj_matrix01, -1).T
                #adj_matrix10 = np.tril(adj_matrix10, -1) + np.tril(adj_matrix10, -1).T

                #print(adj_matrix00)

                adj_matrix = np.zeros((N,N))

                adj_matrix[:gs0,:gs0] = adj_matrix00
                adj_matrix[gs0:,gs0:] = adj_matrix11
                adj_matrix[gs0:,:gs0] = adj_matrix01
                #adj_matrix[:gs0,gs0:] = adj_matrix10
                adj_matrix[:gs0,gs0:] = adj_matrix01.T

                g =  ig.Graph.Adjacency(adj_matrix.tolist(), mode=ig.ADJ_UNDIRECTED)


                if not any(x == 0 for x in g.degree()):
                    success = True
            #print("sum of adj matrix",np.sum(adj_matrix))
            #print("column sum of adj matrix",np.sum(adj_matrix,axis=0))
            #print("row sum of adj matrix",np.sum(adj_matrix,axis=1))
            #vertex_color = [(0,0,0) for i in range(N)]
            #for i in range(N):
            #    if membership[i] == 0:
            #        vertex_color[i] = (1,0,0)
            #    else:
            #        vertex_color[i] = (0,1,0)
            #ig.plot(g,target="sample_network.png",vertex_color=vertex_color)
            g.vs["membership"] = membership
            #print(g.vs["membership"])
        elif (P_layer["type"] == "WS_MS"):
            g = ig.Graph.Watts_Strogatz(dim=P_layer["D"], size=P_layer["L"], nei=P_layer["n_neighbors"], p=0)
            g0 = ig.Graph.Watts_Strogatz(dim=P_layer["D"], size=P_layer["L"], nei=P_layer["n_neighbors"], p=0)

            if P_layer["n_neighbors"] == 1 or P_layer["D"] == 1:
                g = MS_rewire(g,g0,dim = P_layer["D"], p_rewire = P_layer.get("p_rewire",0))

        else:
            print("GRAPH: " + P_layer["type"] + " not implemented yet")

        P_net["N_nodes"] = g.vcount()
        G.append(g)
    return G