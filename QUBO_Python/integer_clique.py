from utils import graph_neighboors, random_edges, draw_graph

def integer_variables(edges: list[tuple], n: int):
    nodes = [i+1 for i in range(n)]
    neighboors = graph_neighboors(edges, n)
    offset = n*'0'
    variables = []
    restrictions = [[] for _ in range(n)]
    for p in range(1,2**n):
        exist = True
        binary = (offset + format(p, 'b'))[-n:]
        new_var = []
        for i in range(n):
            if binary[i-n] == "1":
                new_var.insert(0,nodes[n-1-i])
        if len(new_var) > 1:
            for line in range(len(new_var)-1):
                for col in range(line+1, len(new_var)):
                    if neighboors[new_var[line]-1][new_var[col]-1] != 1:
                        exist = False
        if exist:
            for i in range(n):
                if binary[i-n] == "1":
                    restrictions[nodes[n-1-i]-1].append(p)
            variables.append(new_var)
    return variables, restrictions


if __name__ == "__main__":
    nodes = 23
    edges = random_edges(nodes, 30)
    var, res = integer_variables(edges, nodes)
    print(sorted(var, key= lambda x: len(x)))
    print()
    print(res)
    draw_graph(edges, nodes)