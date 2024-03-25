include("./utils.jl")
using JuMP
using QUBODrivers
using BenchmarkTools
using Graphs
using .RandomGraph
using .Utils
include("./coloring.jl")

n_colors = 2
nodes = 5
edges = [(1,2), (1,3), (1,4), (1,5), (2,3), (2,5), (3,4), (4,5)]

line_graph_nodes = length(edges)
line_graph_edges = line_graph(edges)

run(line_graph_nodes, line_graph_edges, n_colors)
