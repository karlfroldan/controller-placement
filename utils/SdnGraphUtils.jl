module SdnGraphUtils

using Graphs
using GLMakie
using GraphMakie

using SimpleWeightedGraphs

export SDNGraph, node_colors, attack_graph, plot_sdn, simple_network, plot_sdn!

mutable struct SDNGraph
    g :: Union{SimpleGraph, SimpleWeightedGraph}
    node_positions :: AbstractVector{Point2f}
    graph_labels :: AbstractVector{AbstractString}
    attacked_nodes :: AbstractVector{Integer}
    controller_nodes :: AbstractVector{Integer}
end

function SDNGraph(g, node_positions, graph_labels)
    SDNGraph(g, node_positions, graph_labels, [], [])
end

function node_colors(graph :: SDNGraph)
    # Node colors:
	# - blue : controller
	# - red  : attacked
	# - purple : attacked controller
	# - gray : otherwise

	controllers = graph.controller_nodes
	attacked = graph.attacked_nodes
	g = graph.g

	attacked_controllers = intersect(controllers, attacked)
	colored_nodes = union(controllers, attacked)
	non_colored_nodes = setdiff(1:nv(g), colored_nodes)

	controllers = setdiff(controllers, attacked_controllers)
	attacked = setdiff(attacked, attacked_controllers)

	# Color them
	color_it(set, color) = map(v -> (v, color), set)
	
	non_colored_nodes = color_it(non_colored_nodes, colorant"#474749")
	controllers = color_it(controllers, colorant"#3fb3fc")
	attacked = color_it(attacked, colorant"#fc6a7f")
	attacked_controllers = color_it(attacked_controllers, colorant"#b886f9")

	all_nodes = union(non_colored_nodes, controllers, attacked, attacked_controllers)
	sort!(all_nodes; by=x -> x[1])
	map(x -> x[2], all_nodes)
end

function my_add_edge!(g, s, d; w = 0)
    if w != 0
        add_edge!(g, s, d, w)
    else
        add_edge!(g, s, d)
    end
end

function simple_network(weighted = false)
	g = weighted ? SimpleWeightedGraph() : SimpleGraph()

    # Squeeze weights between [100, 200]
    w(i, j) = weighted ? (100 + (100 * log(1 + abs(i - j))) / log(9)) : 0

	add_vertices!(g, 9)

	my_add_edge!(g, 1, 2; w = w(1, 2))
	my_add_edge!(g, 2, 3; w = w(2, 3))
	my_add_edge!(g, 3, 4; w = w(3, 4))
	my_add_edge!(g, 1, 3; w = w(1, 3))
	my_add_edge!(g, 3, 5; w = w(3, 5))
	my_add_edge!(g, 4, 5; w = w(4, 5))
	my_add_edge!(g, 5, 8; w = w(5, 8))
	my_add_edge!(g, 5, 7; w = w(5, 7))
	my_add_edge!(g, 6, 7; w = w(6, 7))
	my_add_edge!(g, 7, 8; w = w(7, 8))
	my_add_edge!(g, 7, 9; w = w(7, 9))

	positions = Point2f[
		(-2.0, 0.0),
		(-1.0, 1.0),
		(-1.0, 0.0),
		(-1.0, -1.0),
		(0.0, 0.0),
		(1.0, 1.0),
		(1.0, 0.0),
		(1.0, -1.0),
		(2.0, 0.0),
	]

	labels = map(x -> "v" * string(x), 1:nv(g))
	
	SDNGraph(g, positions, labels)
end

function plot_sdn(graph :: SDNGraph)
	colors = node_colors(graph)

    offsets = [Point2f(0.1, 0.1) for _ in 1:nv(graph.g)]

	graphplot(
		graph.g; 
		layout=graph.node_positions, 
		nlabels=graph.graph_labels,
        nlabels_offset=offsets,
		node_color=colors,
        edge_color=colorant"#7f7f7f"
	)
end

function plot_sdn!(ax, graph :: SDNGraph)
	colors = node_colors(graph)

    offsets = [Point2f(0.1, 0.1) for _ in 1:nv(graph.g)]

	graphplot!(
        ax,
		graph.g; 
		layout=graph.node_positions, 
		nlabels=graph.graph_labels,
        nlabels_offset=offsets,
		node_color=colors,
        edge_color=colorant"#7f7f7f"
	)
end

end
