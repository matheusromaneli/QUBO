include("./utils.jl")
using JuMP
using QUBODrivers
using BenchmarkTools
using .Utils

p = 2

function r1(nodes:: Int64 , edges:: Vector{Vector{Int64}} , n_tess:: Int64, possible_edges::Vector{Tuple{Int64, Int64}},possible_ancilla::Vector{Tuple{Int64, Int64, Int64}}, Q::Vector{Vector{Int64}})
    edges_vars = length(possible_edges)
    for t_index =1:n_tess
        init = (t_index-1) * edges_vars
        ancilla_init = n_tess*edges_vars
        ancilla_vars = length(possible_ancilla)
        a_init = ancilla_init + (t_index-1)*ancilla_vars*3
        for a_index =1:ancilla_vars
            # Apply restriction: 
            #  2*a*b + 2*a*c + 2*b*c 
            # -2*a*x - 2*b*x + 3*x 
            # -2*a*y - 2*c*y + 3*y 
            # -2*b*z - 2*c*z + 3*z 
            # -a*z - b*y - c*x
            (a,b,c) = possible_ancilla[a_index]
            #Index of my ancilla bits
            x_index = a_init + (a_index-1)*3 + 1
            y_index = a_init + (a_index-1)*3 + 2
            z_index = a_init + (a_index-1)*3 + 3
            #Index of my edges based on a, b, c variables
            pos_e1 = init+findfirst(item -> item == (a,b),possible_edges)
            pos_e2 = init+findfirst(item -> item == (a,c),possible_edges)
            pos_e3 = init+findfirst(item -> item == (b,c),possible_edges)
            # Escolher par com inexistente
            Q[pos_e1][pos_e2] += 2*p
            Q[pos_e1][pos_e3] += 2*p
            Q[pos_e2][pos_e3] += 2*p
            # Penalidade auxiliar para a variável a mais x
            Q[pos_e1][x_index] += -2*p
            Q[pos_e2][x_index] += -2*p
            Q[x_index][x_index] = 3*p
            # Penalidade auxiliar para a variável a mais y
            Q[pos_e1][y_index] += -2*p
            Q[pos_e3][y_index] += -2*p
            Q[y_index][y_index] = 3*p
            # Penalidade auxiliar para a variável a mais z
            Q[pos_e2][z_index] += -2*p
            Q[pos_e3][z_index] += -2*p
            Q[z_index][z_index] = 3*p
            #Penalidade por escolher todas as arestas
            Q[pos_e1][z_index] += -p
            Q[pos_e2][y_index] += -p
            Q[pos_e3][x_index] += -p
            #Auxiliar da seleção de todas arestas
            #TODO: testar como se fosse C(A(1-B)) para não haver termos unicos 
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
                    Q[init+index][init+index] += 2*p
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
    
    size = n_tess * (edges_vars + ancilla_vars*3)

    Q = [[0 for _ = 1:size] for _ = 1:size]
    _r2 = r2(nodes, edges, n_tess, possible_edges,possible_ancilla, Q)
    _r1 = r1(nodes, edges, n_tess, possible_edges,possible_ancilla, Q)
    return reduce(hcat, Q)',size
end

# nodes = 3
# edges = [
#     [0,0,0],
#     [0,0,0],
#     [0,0,0],
# ]
nodes = 4
# square with diagonal (1,3) and not (2,4)
edges = [
    [0,1,1,1],
    [1,0,1,1],
    [1,1,0,1],
    [1,1,1,0],
]
# square
# edges = [
#     [0,1,0,1],
#     [1,0,1,0],
#     [0,1,0,1],
#     [1,0,1,0],
# ]
n_tess = 1

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