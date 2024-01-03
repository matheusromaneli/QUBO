include("./utils.jl")
using JuMP
using QUBODrivers
using BenchmarkTools
using .Utils

p = 4

function r1(nodes:: Int64 , edges:: Vector{Vector{Int64}} , n_tess:: Int64, possible_edges::Vector{Tuple{Int64, Int64}},possible_ancilla::Vector{Tuple{Int64, Int64, Int64}}, Q::Vector{Vector{Int64}})
    edges_vars = length(possible_edges)
    for t_index =1:n_tess
        init = (t_index-1) * edges_vars
        for a1_index =1:edges_vars
            for a2_index =a1_index+1:edges_vars
                a1 = possible_edges[a1_index]
                a2 = possible_edges[a2_index]
                if has_common(a1,a2)
                    #Penalidade por escolher 2 arestas em comum
                    Q[init+a1_index][init+a2_index] = p
                    Q[init+a2_index][init+a1_index] = p
                else
                    Q[init+a1_index][init+a2_index] = -p
                    Q[init+a2_index][init+a1_index] = -p
                end
            end
        end
        ancilla_init = n_tess*edges_vars
        ancilla_vars = length(possible_ancilla)
        a_init = ancilla_init + (t_index-1)*ancilla_vars
        for a_index =1:ancilla_vars
            (x,y,z) = possible_ancilla[a_index]
            pos_e1 = findfirst(item -> item == (x,y),possible_edges)
            pos_e2 = findfirst(item -> item == (x,z),possible_edges)
            pos_e3 = findfirst(item -> item == (y,z),possible_edges)
            e1_signal = (edges[x][y] == 1) ? 1 : -1
            e2_signal = (edges[x][z] == 1) ? 1 : -1
            e3_signal = (edges[y][z] == 1) ? 1 : -1
            #Penalidade por escolher todos os 3 vértices
            Q[a_init+a_index][a_init+a_index] = p
            #Penalidade para cada aresta da tripla
            Q[a_init+a_index][init+pos_e2] = p * e1_signal
            Q[a_init+a_index][init+pos_e1] = p * e2_signal
            Q[a_init+a_index][init+pos_e3] = p * e3_signal
            Q[init+pos_e1][a_init+a_index] = p * e1_signal
            Q[init+pos_e2][a_init+a_index] = p * e2_signal
            Q[init+pos_e3][a_init+a_index] = p * e3_signal
        end
    end
end

function r2(nodes:: Int64 , edges:: Vector{Vector{Int64}} , n_tess:: Int64, possible_edges::Vector{Tuple{Int64, Int64}},possible_ancilla::Vector{Tuple{Int64, Int64, Int64}}, Q::Vector{Vector{Int64}})
    edges_vars = length(possible_edges)
    for t_index = 1:n_tess
        init = (t_index-1) * edges_vars
        for line = 1:nodes
            for column = line+1:nodes
                index = findfirst(item -> item == (line,column),possible_edges)
                if edges[line][column] == 0
                    # Penalidade por escolher não existente
                    Q[init+index][init+index] = p
                else
                    Q[init+index][init+index] = -p
                    for t_aux = 1:n_tess
                        if t_aux != t_index
                            t_init = (t_aux-1) * edges_vars
                            #Penalidade por escolher existente em alguma tesselação
                            Q[init+index][t_init+index] = p/2
                        end
                    end
                end
            end
        end
    end
end
            
function Q(nodes:: Int64 , edges:: Vector{Vector{Int64}} , n_tess:: Int64)
    possible_edges = [(x,y) for x =1:nodes for y =x+1:nodes]
    edges_vars = length(possible_edges)
    possible_ancilla = [(x,y,z) for x =1:nodes for y =1x+1:nodes for z = y+1:nodes]
    ancilla_vars = length(possible_ancilla)
    
    size = n_tess * (edges_vars + ancilla_vars)

    Q = [[0 for _ = 1:size] for _ = 1:size]
    _r1 = r1(nodes, edges, n_tess, possible_edges,possible_ancilla, Q)
    _r2 = r2(nodes, edges, n_tess, possible_edges,possible_ancilla, Q)
    return reduce(hcat, Q)',size
end

nodes = 4
edges = [
    [0,1,1,1],
    [1,0,1,0],
    [0,1,0,1],
    [1,0,1,0],
]
n_tess = 2

q, vars = Q(nodes,edges,n_tess)

show(stdout, "text/plain", q)

@time begin
    model = Model(ExactSampler.Optimizer)

    @variable(model, x[1:(vars)], Bin)
    @objective(model, Min, x' * q * x)

    optimize!(model)
end

xi = value.(x)
show(xi)