using JuMP
using QUBODrivers

p = 1

function r1(nodes:: Int64, edges:: Vector{Tuple{Int64, Int64}}, Q::Vector{Vector{Int64}})
    for (neigh1, neigh2) in edges
        Q[neigh1][neigh1] += -p
        Q[neigh2][neigh2] += -p
    end
    show(stdout, "text/plain", Q)
end

function r2(nodes:: Int64, edges:: Vector{Tuple{Int64, Int64}}, Q::Vector{Vector{Int64}})
    triples = [(x,y,z) for x =1:nodes for y = x+1:nodes for z = y+1:nodes]
    for (a,b,c) in triples
        pos_e1 = findfirst(item -> item == (a,b), edges)
        pos_e2 = findfirst(item -> item == (a,c), edges)
        pos_e3 = findfirst(item -> item == (b,c), edges)
        
        if pos_e1 !== nothing && pos_e2 !== nothing && pos_e3 !== nothing
            println("$(a) $(b) $(c)")
            Q[a][b] += p
            Q[a][c] += p
            Q[b][c] += p
        end
    end
end

function Q(nodes:: Int64 , edges:: Vector{Tuple{Int64, Int64}})
    Q = [[0 for _ = 1:nodes] for _ = 1:nodes]
    r1(nodes, edges, Q)
    r2(nodes, edges, Q)

    return reduce(hcat, Q)', nodes
end

nodes = 6
edges = [(1,2), (1,3), (1,4), (1,5), (1,6), (2,3), (3,4)]

q, vars = Q(nodes, edges)

show(stdout, "text/plain", q)

model = Model(ExactSampler.Optimizer)

@variable(model, x[1:nodes], Bin)
@objective(model, Min, (x' * q * x))

optimize!(model)

xi = value.(x; result=1)
yi = objective_value(model; result=1)
println(xi)
println("Value:", yi)