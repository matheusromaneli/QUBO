include("./utils.jl")
include("./coloring.jl")
using JuMP
using QUBODrivers
using BenchmarkTools
using Graphs
using .RandomGraph
using .Utils

n_colors = 4

nodes = 5
edges = [(1,2), (1,3), (1,4), (1,5), (2,3), (2,5), (3,4), (4,5)]

graph_edges = Edge.(edges)
graph = SimpleGraph(graph_edges)
clique_graph = maximal_cliques(graph)

clique_graph_nodes = length(clique_graph)
clique_graph_edges_result = clique_graph_edges(clique_graph, clique_graph_nodes)

run(clique_graph_nodes, clique_graph_edges_result, n_colors)