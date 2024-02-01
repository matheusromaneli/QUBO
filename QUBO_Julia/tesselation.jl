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
        ancilla_init = n_tess*edges_vars
        ancilla_vars = length(possible_ancilla)
        a_init = ancilla_init + (t_index-1)*ancilla_vars
        for a_index =1:ancilla_vars
            (x,y,z) = possible_ancilla[a_index]
            pos_e1 = findfirst(item -> item == (x,y),possible_edges)
            pos_e2 = findfirst(item -> item == (x,z),possible_edges)
            pos_e3 = findfirst(item -> item == (y,z),possible_edges)
            # Penalidade auxiliar para a variável a mais
            Q[init+pos_e1][init+pos_e2] += p/2
            Q[init+pos_e1][a_init+a_index] += -2*p
            Q[init+pos_e2][a_init+a_index] += -2*p
            Q[a_init+a_index][a_init+a_index] = 3*p
            # Escolher 2 arestas apenas
            Q[init+pos_e1][init+pos_e2] += p/2
            Q[init+pos_e1][init+pos_e3] += p/2
            Q[init+pos_e2][init+pos_e3] += p/2
            # Escolher par com inexistente
            if edges[x][y] == 0
                Q[init+pos_e1][init+pos_e2] += p
                Q[init+pos_e1][init+pos_e3] += p
            end
            if edges[x][z] == 0
                Q[init+pos_e1][init+pos_e2] += p
                Q[init+pos_e2][init+pos_e3] += p
            end
            if edges[y][z] == 0
                Q[init+pos_e2][init+pos_e3] += p
                Q[init+pos_e1][init+pos_e3] += p
            end
            #Penalidade por escolher todas as arestas
            Q[init+pos_e3][a_init+a_index] += -3*p
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
                    Q[init+index][init+index] += p
                else
                    # Penalidade por escolher existente
                    Q[init+index][init+index] += -p
                    for t_aux = t_index+1:n_tess
                        if t_aux != t_index
                            t_init = (t_aux-1) * edges_vars
                            #Penalidade por escolher existente em alguma tesselação
                            Q[init+index][t_init+index] += p/2
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
    _r2 = r2(nodes, edges, n_tess, possible_edges,possible_ancilla, Q)
    _r1 = r1(nodes, edges, n_tess, possible_edges,possible_ancilla, Q)
    return reduce(hcat, Q)',size
end

# nodes = 3
# #square with diagonal (1,3) and not (2,4)
# edges = [
#     [0,1,0],
#     [1,0,0],
#     [0,0,0],
# ]
nodes = 4
#square with diagonal (1,3) and not (2,4)
edges = [
    [0,1,1,1],
    [1,0,1,0],
    [1,1,0,1],
    [1,0,1,0],
]
# square
# edges = [
#     [0,1,0,1],
#     [1,0,1,0],
#     [0,1,0,1],
#     [1,0,1,0],
# ]
n_tess = 2

q, vars = Q(nodes,edges,n_tess)

show(stdout, "text/plain", q)

# q = [
#      [1,-1,0,-1],
#      [0,0,0,1],
#      [0,0,1,-1],
#      [0,0,0,0],
#     ]
# q = reduce(hcat, q)'
# vars=4

@time begin
    model = Model(ExactSampler.Optimizer)

    @variable(model, x[1:(vars)], Bin)
    @objective(model, Min, 6*p + x' * q * x)

    optimize!(model)
end

println()
possible_edges = [(x,y) for x =1:nodes for y =x+1:nodes]
edges_vars = length(possible_edges)
show(possible_edges)
println()
for i = 1:2
    xi = value.(x; result=i)
    yi = objective_value(model; result=i)
    println(xi)
    println("Value:", yi)
    for i=1:(edges_vars*n_tess)
        if xi[i] == 1.0
            println("Edge ", possible_edges[i - (edges_vars*div(i-1,edges_vars))], " in tesselation ", 1+div(i,(edges_vars+1)))
        end
    end
    println()
end