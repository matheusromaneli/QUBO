include("coloring.jl")
using Graphs

nodes = 5
all_edges = [(1,2), (1,3), (1,4), (1,5), (2,3), (2,5), (3,4), (4,5)]

biggest_coloring = 0
for v = 1:nodes
    # N(v)
    # arestas diretamente ligada a v
    neighbors_nodes_v::Vector{Integer} = []
    for (v1,v2) in all_edges
        if v==v1
            push!(neighbors_nodes_v, v2)
        elseif v==v2
            push!(neighbors_nodes_v, v1)
        end
    end
    graph_edges = Edge.(all_edges)
    graph = SimpleGraph(graph_edges)

    neighbors, vmap = induced_subgraph(graph, neighbors_nodes_v)

    final_graph = complement(neighbors)

    translate_graph:: Vector{Tuple{Int64, Int64}}= []

    foreach(item-> push!(translate_graph, (vmap[item.src], vmap[item.dst])), collect(edges(final_graph)))

    successfull = false
    edge_colors = []
    i=1
    while successfull == false
        if i < 4
            i += 1
            edge_colors, successfull = run(nodes, translate_graph, i)
        else
            successfull = true
            edge_colors = "Need more than 4"
        end
    end
    println("complement = $translate_graph")
    println("Node $v is successfull with $i colors:\n $edge_colors\n")
    if i > biggest_coloring
        global biggest_coloring = i
    end
end

println("Biggest value is: $biggest_coloring")