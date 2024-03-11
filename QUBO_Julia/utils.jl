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
        return edges
    end
    export random_edges
end

module Utils
    function node_in_common(nodes_1:: Vector{Int64}, nodes_2:: Vector{Int64})
        for n in nodes_1
            for m in nodes_2
                if n == m
                    return true
                end
            end
        end
        return false
    end

    function clique_graph_edges(clique_graph:: Vector{Vector{Int64}}, nodes:: Int64)
        edges::Vector{Tuple{Int64,Int64}} = []
        for i = 1:nodes
            for j = i+1:nodes
                if node_in_common(clique_graph[i], clique_graph[j])
                    push!(edges, (i,j))
                end
            end
        end
        return edges
    end

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
    export has_common, vector_to_matrix, clique_graph_edges
end
