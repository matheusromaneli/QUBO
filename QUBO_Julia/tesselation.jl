include("./utils.jl")
using JuMP
using QUBODrivers
using BenchmarkTools
using Graphs
using .Utils

p = 2

function certificate(n_tess::Int64, nodes:: Int64, adjacency_matrix:: Vector{Vector{Int64}}, choosen:: Vector{Float64}):: Bool
    edge_list = adjacency_to_list(adjacency_matrix)
    edge_choosen = [0 for _ = 1:length(edge_list)]
    possible_edges = [(x,y) for x =1:nodes for y =x+1:nodes]
    #-Verificar cobertura de arestas
    #para cada aresta
    for e in eachindex(edge_list)
        #-Verificar se foi escolhida em alguma tesselação
        #indice na lista de arestas possiveis
        edge_index = findfirst(item -> item == edge_list[e], possible_edges)
        #para cada tesselação
        for t = 1:n_tess
            #adiciono um contador em uma lista auxiliar para saber se foi escolhida
            result_index = edge_index + (t-1) * length(possible_edges)
            if choosen[result_index] == 1.0
                edge_choosen[e] = 1
            end
        end
    end
    #-Verificar lista auxiliar
    #para cada elemento
    for element in edge_choosen
        #se não for escolhido retornar false
        if element == 0
            # println("wrong because doesnt cover all edges")
            return false
        end
    end
    #-Verificar cobertura por tesselações
    #para cada tesselação
    for t = 1:n_tess
        #gerar lista de arestas
        start = (t-1)*length(possible_edges) + 1
        finish= (t) * length(possible_edges)
        tess_choosen = choosen[start:finish]
        tess_edges = []
        for i in eachindex(tess_choosen)
            if tess_choosen[i] == 1.0
                push!(tess_edges, possible_edges[i])
            end
        end
        #-Verificar tesselação
        #para cada aresta da tesselação
        for pivot_index in eachindex(tess_edges)
            pivot = tess_edges[pivot_index]
            #procurar por adjacente
            for aux_index = pivot_index+1:length(tess_edges)
                aux = tess_edges[aux_index]
                #se houver adjacente
                if has_common(pivot, aux)
                    #procurar por complemento
                    completer = third_edge(pivot, aux)
                    found = false
                    #se não houver complemento
                    for edge in tess_edges
                        if completer == edge
                            found = true
                        end
                    end
                    # retorna falso
                    if found == false
                        # println("wrong because have adjacent vertex without a clique")
                        return false
                    end
                end
            end
        end
    end
    #retorna verdadeiro        
    return true
end

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
                end
            end
        end
    end
end
         
function r3(nodes:: Int64 , edges:: Vector{Vector{Int64}} , n_tess:: Int64, possible_edges::Vector{Tuple{Int64, Int64}},possible_ancilla::Vector{Tuple{Int64, Int64, Int64}}, Q::Vector{Vector{Int64}})
    edges_vars = length(possible_edges)
    ancilla_vars = length(possible_ancilla)
    for line = 1:nodes
        for column = line+1:nodes
            index = findfirst(item -> item == (line,column),possible_edges)
            for t_index_i = 1:n_tess
                init_i = (t_index_i-1) * edges_vars
                if edges[line][column] == 1
                    # Penalidade por escolher existente
                    Q[init_i+index][init_i+index] += -p
                    for t_index_j = t_index_i+1:n_tess
                        init_j = (t_index_j-1) * edges_vars
                        #Penalidade por escolher existente em duas tesselação
                        Q[init_i+index][init_j+index] += p
                        init_r3_ancilla = n_tess * (edges_vars + ancilla_vars*3)
                        for t_index_k = t_index_j+1:n_tess
                            init_k = (t_index_k-1) * edges_vars
                            Q[init_i+index][init_j+index] += p
                            Q[init_i+index][init_r3_ancilla+index] -= 2*p
                            Q[init_j+index][init_r3_ancilla+index] -= 2*p
                            Q[init_r3_ancilla+index][init_r3_ancilla+index] += 3*p
                            Q[init_k+index][init_r3_ancilla+index] -=p
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
    r2_ancilla = [(x,y,z) for x =1:nodes for y = x+1:nodes for z = y+1:nodes]
    r2_ancilla_vars = length(r2_ancilla)
    
    r3_ancilla = [(x,y,z) for x =1:n_tess for y = x+1:n_tess for z = y+1:n_tess]
    r3_ancilla_vars = length(r3_ancilla)

    size = n_tess * (edges_vars + r2_ancilla_vars*3) + r3_ancilla_vars*edges_vars

    Q = [[0 for _ = 1:size] for _ = 1:size]
    _r1 = r1(nodes, edges, n_tess, possible_edges, r2_ancilla, Q)
    _r2 = r2(nodes, edges, n_tess, possible_edges, r2_ancilla, Q)
    _r3 = r3(nodes, edges, n_tess, possible_edges, r2_ancilla, Q)

    return reduce(hcat, Q)',size
end

nodes = 3
edges = [
    [0,1,0],
    [1,0,1],
    [0,1,0],
]
# nodes = 4
# # square with diagonal (2,4) and not (1,3)
# edges = [
#     [0,1,0,1],
#     [1,0,1,1],
#     [0,1,0,1],
#     [1,1,1,0],
# ]
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

ancilla_bits = [(x,y,z) for x =1:nodes for y =1x+1:nodes for z = y+1:nodes]
println()
show(ancilla_bits)
println()
println(vars)
# exit(1)

model = Model(ExactSampler.Optimizer)

@variable(model, x[1:(vars)], Bin)
@objective(model, Min, (x' * q * x))

optimize!(model)

println()
possible_edges = [(x,y) for x =1:nodes for y =x+1:nodes]
ancilla_bits = [(x,y,z) for x =1:nodes for y =1x+1:nodes for z = y+1:nodes]
edges_vars = length(possible_edges)
show(possible_edges)
println()
show(ancilla_bits)
println("\n")

if certificate(n_tess, nodes, edges, value.(x; result=1))
    println("Successfull solution!")
    xi = value.(x; result=1)
    yi = objective_value(model; result=1)
    println(xi)
    println("Value:", yi)
    for i=1:(edges_vars*n_tess)
        if xi[i] == 1.0
            println("Edge ", possible_edges[i - (edges_vars*div(i-1,edges_vars))], " in tesselation ", 1+div(i,(edges_vars+1)))
        end
    end
    println()
else
    println("Unsuccessfull solution!")
    xi = value.(x; result=1)
    yi = objective_value(model; result=1)
    println(xi)
    println("Value:", yi)
    for i=1:(edges_vars*n_tess)
        if xi[i] == 1.0
            println("Edge ", possible_edges[i - (edges_vars*div(i-1,edges_vars))], " in tesselation ", 1+div(i,(edges_vars+1)))
        end
    end
    println()
end