import numpy as np
import igraph as ig


def n_simple_rewire(N,p, size = 1):
    # calculates what to put for   n   in g.rewire(n = ?, mode = "simple")   as to be consistent with the Watts-Strogatz algorithm
    # this function affects exactly and determinstically    2*n    edges.
    # In the Watts Strogatz model,  for   N    edges    and a rewiring probability   p   , on average p*N edges will be affected.
    # The exact value comes from the probability distribution    binomial(N,p)
    # The following formula captures the stochastic nature and is close to half the average of the above distribution.
    # close means: off by less than 0.25  (so the total number of affected edges is off by less than 0.5)

    #N:     number of edges
    #p:     rewiring probability
    #size:  how many values to generate

    return (    0.5* np.random.binomial ( n = N, p=p, size = size )            ).astype(int)

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
        P_layer = P_net["Layer_" + str(i)]
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

            #Code to display all neighbors for each node, to check that the adj matrix is working

        elif(P_layer["type"] == "Moore_Lattice_WS_MS"):
            #a lattice with Moore neighborhood, i.e. 8 neighbors, then rewired Maslov-Sneppen-like
            #a Cellular automaton
            #parameters are "Lx", "Ly", "circular"

            A = MooreLattice_adjacency(P_layer["Lx"],P_layer["Ly"],P_layer.get("circular",False))
            g = ig.Graph.Adjacency(A, mode="undirected")
            g.rewire(n = n_simple_rewire(N = g.ecount(), p=P_layer["P"]) , mode ="simple")
            #g.rewire(n = 1):   two edges are picked and crosswired.
            # e.g.  we pick the edges    0--312; 547--548
            # we delete those edges and add either     0--547; 312--548     or    0--548; 312--547
            # consistency with Watts-Strogatz definition:  there for 200 edges and p=5%  we want to rewire 10 edges.  that means we set n = 0.5 * 10 = 5
            
        elif(P_layer["type"] == "WS"):
            if(np.power(P_layer["L"],P_layer["D"]) != N):
                print("GRAPH SIZE WARNING: L^D != N: " + str(np.power(P_layer["L"],P_layer["D"])) + " != " + str(N))
            g = ig.Graph.Watts_Strogatz(dim=P_layer["D"], size=P_layer["L"], nei=P_layer["NFN"], p=P_layer["P"])
        elif P_layer["type"] == "WS_MS":
            #Watts-Strogatz network with Maslov-Sneppen-like rewiring
            if(np.power(P_layer["L"],P_layer["D"]) != N):
                print("GRAPH SIZE WARNING: L^D != N: " + str(np.power(P_layer["L"],P_layer["D"])) + " != " + str(N))
            g = ig.Graph.Watts_Strogatz(dim=P_layer["D"], size=P_layer["L"], nei=P_layer["NFN"], p=0)
            g.rewire(n = n_simple_rewire(N = g.ecount(), p=P_layer["P"]) , mode ="simple")
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
        else:
            print("GRAPH: " + P_layer["type"] + " not implemented yet")

        P_net["N_nodes"] = g.vcount()
        G.append(g)
    return G