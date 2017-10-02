# Times network inference for different numbers of genes (nodes)
function get_times_per_number_of_genes(algorithm, discretizer, number_of_cells, min_number_of_nodes,
        max_number_of_nodes, step)

    println("Generating data...")
    data_dictionary = Dict()
    for i in 1 : max_number_of_nodes
        binned_values = zeros(Int, number_of_cells)
        number_of_bins = get_bin_ids!(randn(number_of_cells), discretizer, Int(round(sqrt(number_of_cells))), binned_values)
        get_frequencies_from_bin_ids(binned_values, number_of_bins)
        probabilities = get_probabilities("maximum_likelihood", get_frequencies_from_bin_ids(binned_values, number_of_bins))
        data_dictionary[string(i)] = (binned_values, number_of_bins, probabilities)
    end

    println("Storing data as Nodes..")
    nodes = Array{Node}(max_number_of_nodes)
    for (i, x) in enumerate(data_dictionary)
        nodes[i] = Node(x[1], x[2][1], x[2][2], x[2][3])
    end

    println("Timing network inference...")
    sizes = []
    times = []
    for i in min_number_of_nodes : step : max_number_of_nodes
        println("-------------")
        println("$i genes:")
        time = @elapsed InferredNetwork(algorithm, nodes[1:i])
        println("$time s")
        push!(sizes, i)
        push!(times, time)
    end
    
    return sizes, times

end

# Times network inference for different numbers of cells
function get_times_per_number_of_cells(algorithm, discretizer, number_of_nodes, min_number_of_cells,
    max_number_of_cells, step)

    # Infer the network
    function get_network(i)
        nodes = Array{Node}(number_of_nodes)
        for (j, x) in enumerate(data_dictionary)
            binned_values = zeros(Int, i)
            data = x[2][1:i]
            number_of_bins = get_bin_ids!(data, discretizer, 0, binned_values)
            probabilities = get_probabilities("maximum_likelihood", get_frequencies_from_bin_ids(binned_values, number_of_bins))
            nodes[j] = Node(x[1], binned_values, number_of_bins, probabilities)
        end
        InferredNetwork(MINetworkInference(), nodes)
    end

    println("Generating data...")
    data_dictionary = Dict()
    for i in 1 : number_of_nodes
        data_dictionary[string(i)] = randn(max_number_of_cells)
    end


    println("Timing discretization and network inference...")
    sizes = []
    times = []
    for i in min_number_of_cells : step : max_number_of_cells
        println("-------------")
        println("Cells: $i")
        time = @elapsed get_network(i)
        println("$time s")
        push!(sizes, i)
        push!(times, time)
    end

    return sizes, times

end
