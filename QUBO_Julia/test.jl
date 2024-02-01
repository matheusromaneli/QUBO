using JuMP
using QUBODrivers
using BenchmarkTools

vars = 8
p=2
Q = [
    [4,8 , 6,2,0 , 0,-8, 0],
    [0  ,-4, 2,0,2 , 0,-8, 0],
    [0,0 , -4,0,0 , 0 , -12, 0],
    [0,0 , 0 , 4,8 , 6, 0,-8],
    [0,0 , 0,0 , -4,  2, 0,-8],
    [0,0 , 0,0,0 , -4, 0 , -12],
    [0,0 , 0,0,0 , 0,12, 0],
    [0,0 , 0,0,0 , 0, 0,12],
]

q = reduce(hcat, Q)'
show(q)
println()
@time begin
    model = Model(ExactSampler.Optimizer)

    @variable(model, x[1:(vars)], Bin)
    @objective(model, Min, x' * q * x)

    optimize!(model)
end

for i = 1:result_count(model)
    xi = value.(x; result=i)
    yi = objective_value(model; result=i)

    println("f($xi) = $yi")
end