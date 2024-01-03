import networkx as nx
from random import randint
from matplotlib import pyplot as plt

def graph_neighboors(edges: list[tuple], n: int):
    '''Return a matrix which pos[i][j] is 1 if node_i has conection with node_j and vice-versa'''
    adj_list = [[0 for _ in range(n)] for _ in range(n)]
    for (src,dest) in edges:
        adj_list[src-1][dest-1] = 1
        adj_list[dest-1][src-1] = 1
    return adj_list

def has_common(a1: tuple[int,int],a2: tuple[int,int]):
    for a in a1:
        for b in a2:
            if a == b:
                return True
    return False

def draw_graph(edges: list[tuple], n: int, color_map: list[str] = None):
    '''Draw graph'''
    G = nx.Graph()
    G.add_nodes_from([i for i in range(1,n+1)])
    G.add_edges_from(edges)
    nx.draw(G, node_color = color_map, with_labels = True)
    plt.show()

def is_pow_of_two(n: int):
    return (n != 0) & (n & (n-1) == 0)

def random_edges(num_verts: int, edge_ratio: int):
    edges = []
    for i in range(1,num_verts+1):
        for j in range(1,num_verts+1):
            if i!=j and randint(1,100) < edge_ratio: edges.append((i,j)) 
    return edges

def show(matrix: list[list]):
    for line in matrix:
        print(line)
