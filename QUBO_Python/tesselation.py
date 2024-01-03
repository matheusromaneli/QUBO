from utils import has_common

vert = 4
edges = [
    [0,1,0,1],
    [1,0,1,0],
    [0,1,0,1],
    [1,0,1,0],
]
tess = 2
possible_edges = [(x,y) for x in range(vert) for y in range(x+1,vert)]
print(possible_edges)
edges_vars = len(possible_edges)
ancilla_init = tess*edges_vars
possible_ancilla = [(x,y,z) for x in range(vert) for y in range(x+1,vert) for z in range(y+1,vert)]
print(possible_ancilla)
ancilla_vars = len(possible_ancilla)

matrix_size = tess * (edges_vars + ancilla_vars)
a = [[0 for _ in range(matrix_size)] for _ in range(matrix_size)]
def r1(m):
    for t_index in range(tess):
        init = t_index * edges_vars
        for a1_index in range(edges_vars):
            for a2_index in range(a1_index+1, edges_vars):
                a1 = possible_edges[a1_index]
                a2 = possible_edges[a2_index]
                if has_common(a1,a2):
                    m[init+a1_index][init+a2_index] = 1
                    m[init+a2_index][init+a1_index] = 1

        a_init = ancilla_init + t_index*ancilla_vars
        for a_index in range(ancilla_vars): 
            (x,y,z) = possible_ancilla[a_index]
            e1 = possible_edges.index((x,y))
            e2 = possible_edges.index((x,z))
            e3 = possible_edges.index((y,z))
            m[a_init+a_index][a_init+a_index] = -3
            m[a_init+a_index][init+e1] = 1
            m[a_init+a_index][init+e2] = 1
            m[a_init+a_index][init+e3] = 1
            m[init+e1][a_init+a_index] = 1
            m[init+e2][a_init+a_index] = 1
            m[init+e3][a_init+a_index] = 1

def r2and3(m):
    for t_index in range(tess):
        init = t_index * edges_vars
        for line in range(vert):
            for column in range(line+1,vert):
                index = possible_edges.index((line,column))
                if edges[line][column] == 0:
                    m[init+index][init+index] = 1
                else:
                    m[init+index][init+index] = -1
                    for t_aux in range(tess):
                        t_init = t_aux * edges_vars
                        if t_index != t_aux:
                            m[init+index][t_init+index] = 1
            
r1(a)
r2and3(a)
for line in a:
    print(line)