include("./utils.jl")
using JuMP
using QUBODrivers
using BenchmarkTools
using .RandomGraph
using .Utils

function certificate(edges:: Vector{Tuple{Int64,Int64}}, colors:: Vector{Int64})
    for (e1,e2) in edges
        if colors[e1] == colors[e2]
            return false
        end
    end
    return true
end

function r1(nodes:: Int64, edges:: Vector{Tuple{Int64,Int64}} , n_colors:: Int64, Q::Vector{Vector{Int64}})
    for i= 0:nodes-1
        start = i*n_colors
        for j = start+1:start+n_colors
            for k = start+1:start+n_colors
                if j != k 
                    Q[j][k] = 4 
                else 
                    Q[j][k] =-4
                end
            end
        end
    end
end
function r2(nodes:: Int64 , edges:: Vector{Tuple{Int64,Int64}} , n_colors:: Int64, Q::Vector{Vector{Int64}})
    for (i , j) in edges
        for k = 0:n_colors-1
            m = (i) * n_colors - k
            o = (j) * n_colors - k
            Q[m][o] = 2
            Q[o][m] = 2
        end
    end
end

function Q(nodes:: Int64 , edges:: Vector{Tuple{Int64,Int64}} , n_colors:: Int64)
    size = nodes * n_colors
    Q = [[0 for _ = 1:size] for _ = 1:size]
    _r1 = r1(nodes, edges, n_colors, Q)
    _r2 = r2(nodes, edges, n_colors, Q)
    return reduce(hcat, Q)'
end

function run(nodes:: Int64, edges:: Vector{Tuple{Int64,Int64}}, n_colors:: Int64)
    q = Q(nodes, edges, n_colors)
    model = Model(ExactSampler.Optimizer)

    @variable(model, x[1:(nodes*n_colors)], Bin)
    @objective(model, Min, x' * q * x)

    optimize!(model)

    xi = value.(x)
    nodes_colors = [0 for _ = 1:nodes]
    for node=1:nodes
        for color=0:n_colors-1
            node_index = (node-1)*n_colors + 1
            if xi[node_index + color] == 1.0
                nodes_colors[node] = color
            end    
        end
    end

    successfull = false
    if certificate(edges, nodes_colors)
        successfull = true
    end
    return nodes_colors, successfull
end

nodes = 5
all_edges = random_edges(nodes, 60)
n_colors = 3

run(nodes, all_edges, n_colors)

export run