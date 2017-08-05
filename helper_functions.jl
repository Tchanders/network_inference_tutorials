# Makes Junet graph from a NetworkInference network analysis
# Network analyses rank all possible edges, so also need a percentage of edges to keep
function network_analysis_to_graph(network::NetworkAnalysis, percentage_threshold::Float64)

    # Make empty graph
    graph = Graph()
    graph[:, :label] = [""]

    # Add nodes to graph
    addnode!(graph, length(network.genes))
    
    # Make dict of gene names to incremental Int ids
    genes_dict = Dict{String, Int}()
    for (i, gene) in enumerate(network.genes)
        genes_dict[gene.name] = i
        graph[:, :label][i] = gene.name
    end
    
    # For each edge we want to keep, add the edge to the graph
    number_of_edges = Int(round(percentage_threshold * length(network.edges)))
    for edge in network.edges[1:number_of_edges]
        genes = collect(edge.genes)
        addedge!(graph, genes_dict[genes[1].name], genes_dict[genes[2].name])
    end
    
    return graph
end

# Times network inference for different numbers of genes
function get_times_per_number_of_genes(algorithm, discretizer, number_of_cells, min_number_of_genes,
        max_number_of_genes, step)

    println("Generating data...")
    data_dictionary = Dict()
    for i in 1 : max_number_of_genes
        discretized_values = zeros(Int, number_of_cells)
        number_of_bins = get_bin_ids!(randn(number_of_cells), discretizer, Int(round(sqrt(number_of_cells))), discretized_values)
        get_frequencies_from_bin_ids(discretized_values, number_of_bins)
        probabilities = get_probabilities("maximum_likelihood", get_frequencies_from_bin_ids(discretized_values, number_of_bins))
        data_dictionary[string(i)] = (discretized_values, number_of_bins, probabilities)
    end

    println("Storing data as Genes..")
    genes = Array{Gene}(max_number_of_genes)
    for (i, x) in enumerate(data_dictionary)
        genes[i] = Gene(x[1], x[2][1], x[2][2], x[2][3])
    end

    println("Timing network inference...")
    sizes = []
    times = []
    for i in min_number_of_genes : step : max_number_of_genes
        println("-------------")
        println("$i genes:")
        time = @elapsed NetworkAnalysis(algorithm, genes[1:i])
        println("$time s")
        push!(sizes, i)
        push!(times, time)
    end
    
    return sizes, times

end

# Times network inference for different numbers of genes
function get_times_per_number_of_cells(algorithm, discretizer, number_of_genes, min_number_of_cells,
    max_number_of_cells, step)

    # Infer the network
    function get_network(i)
        genes = Array{Gene}(number_of_genes)
        for (j, x) in enumerate(data_dictionary)
            discretized_values = zeros(Int, i)
            data = x[2][1:i]
            number_of_bins = get_bin_ids!(data, discretizer, 0, discretized_values)
            probabilities = get_probabilities("maximum_likelihood", get_frequencies_from_bin_ids(discretized_values, number_of_bins))
            genes[j] = Gene(x[1], discretized_values, number_of_bins, probabilities)
        end
        NetworkAnalysis(MINetworkInference(), genes)
    end

    println("Generating data...")
    data_dictionary = Dict()
    for i in 1 : number_of_genes
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
