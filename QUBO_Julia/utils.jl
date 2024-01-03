module RandomGraph
    using Random

    function random_edges(num_verts:: Int64, edge_ratio:: Int64)
        edges::Vector{Tuple{Int64, Int64}} = []
        for i = 1:num_verts
            for j = 1:num_verts
                if i!=j && rand(1:100) < edge_ratio
                    push!(edges,(i,j)) 
                end
            end
        end
        println(edges)
        return edges
    end
    export random_edges
end

module Utils
    function vector_to_matrix(m)
        new_vec = []
        for line in m
            for col in line
                push!(new_vec,col)    
            end
        end  
        return new_vec  
    end

    function has_common(edge1:: Tuple{Int64, Int64}, edge2:: Tuple{Int64, Int64})
        for i in edge1
            for j in edge2
                if i == j
                    return true
                end
            end
        end
        return false
    end
    export has_common, vector_to_matrix
end
