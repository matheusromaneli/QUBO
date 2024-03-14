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

    function line_graph(edges_list:: Vector{Tuple{Int64, Int64}}):: Vector{Tuple{Int64, Int64}}
        result:: Vector{Tuple{Int64, Int64}} = []
        for pivot_index in eachindex(edges_list)
            pivot = edges_list[pivot_index]
            for aux_index in pivot_index+1:length(edges_list)
                aux = edges_list[aux_index]
                if has_common(pivot, aux)
                    push!(result, (pivot_index,aux_index))
                end
            end
        end
        return result
    end

    function adjacency_to_list(adjacency_matrix:: Vector{Vector{Int64}}):: Vector{Tuple{Int64,Int64}}
        edges = []
        nodes = length(adjacency_matrix)
        for line = 1:nodes
            for column = line:nodes
                if adjacency_matrix[line][column] == 1
                    push!(edges, (line,column))
                end
            end
        end  
        return edges
    end

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

    function third_edge(edge1:: Tuple{Int64, Int64}, edge2:: Tuple{Int64, Int64}):: Tuple{Int64, Int64}
        if edge1[1] == edge2[1]
            return (edge1[2],edge2[2])
        elseif edge1[1] == edge2[2]
            return (edge1[2],edge2[1])
        elseif edge1[2] == edge2[1]
            return (edge1[1],edge2[2])
        elseif edge1[2] == edge2[2]
            return (edge1[1],edge2[1])
        end
        return (-1,-1)
    end
    export has_common, vector_to_matrix, clique_graph_edges, adjacency_to_list, third_edge, line_graph
end
