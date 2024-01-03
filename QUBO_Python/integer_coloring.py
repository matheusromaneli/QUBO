from utils import random_edges, draw_graph, show
from docplex.mp.model import Model
from qiskit.primitives import Sampler, Estimator
from qiskit.algorithms.optimizers import COBYLA
from qiskit_optimization.translators import from_docplex_mp
from qiskit_optimization.algorithms import MinimumEigenOptimizer
from qiskit import Aer
from qiskit.algorithms import QAOA
from qiskit.utils import QuantumInstance,algorithm_globals
from datetime import datetime


def r1(nodes: int , edges: list [ tuple ] , n_colors: int, Q):
    for i in range (nodes) :
        start = i*n_colors
        for j in range(start, start+n_colors):
            for k in range(start, start+n_colors):
                Q[j][k] = 4 if j != k else -4

def r2(nodes: int , edges: list [ tuple ] , n_colors: int, Q):
    constraints = []
    for (i , j ) in edges :
        model = [1 , "<=" ]
        for k in range ( n_colors ) :
            m = (i-1) * n_colors + k
            o = (j-1) * n_colors + k
            Q[m][o] = 2
            Q[o][m] = 2
            model.append (( i , k ) )
            model.append (( j , k ) )
        constraints.append ( model )
    return constraints

def Q(nodes : int , edges : list [ tuple ] , n_colors : int):
    size = nodes * n_colors
    Q = [[0 for _ in range(size)] for _ in range(size)]
    _r1 = r1(nodes, edges, n_colors, Q)
    _r2 = r2(nodes, edges, n_colors, Q)
    show(Q)
    return Q

def assign_color(result, n_colors):
    colors = ['red', 'yellow', 'green', 'voilet', 'blue'] 
    color_map = []
    for pos, val in enumerate(result.x):
        if val == 1:
            color = pos%n_colors
            color_map.append(colors[color])
    return color_map


if __name__ == "__main__":
    nodes = 5
    color = 3
    edges = random_edges(nodes, 50)
    print()
    q = Q(nodes, edges, color)

    mdl = Model(name="gcolor")
    x = [mdl.binary_var('x%s' % i) for i in range(len(q))]
    objective = mdl.sum([x[i]*q[i][j]* x[j] for i in range(len(q)) for j in range(len(q))])
    mdl.minimize(objective)
    qp = from_docplex_mp(mdl)
    print( qp.prettyprint())
    
    seed = 1234
    t_init = datetime.now()
    algorithm_globals.random_seed = seed
    # sampler = Sampler()
    # optimizer = COBYLA()
    # qaoa = QAOA(sampler, optimizer, reps=2)
    # result = qaoa.compute_minimum_eigenvalue(q)
    quantum_instance = QuantumInstance(backend=Aer.get_backend('qasm_simulator'), shots=1000)
    qaoa = MinimumEigenOptimizer(min_eigen_solver=QAOA(reps=3, quantum_instance=quantum_instance))
    qaoa_result = qaoa.solve(qp)
    t_end = datetime.now()
    print(qaoa_result)
    # print(qaoa_result)
    # color_map = assign_color(qaoa_result, color)
    # while(len(color_map) < nodes):
    #     color_map.append('black')
    # draw_graph(edges, nodes, color_map)

    final = t_end - t_init
    print("Ttotal: ", final)
