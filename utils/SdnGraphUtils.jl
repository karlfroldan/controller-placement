module SdnGraphUtils

using Graphs
using GLMakie
using GraphMakie

using SimpleWeightedGraphs

export SDNGraph, node_colors, attack_graph, plot_sdn, simple_network, plot_sdn!

mutable struct SDNGraph
    g :: Graphs.SimpleGraphs.AbstractSimpleGraph
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

function simple_network(weighted = false)
	g = weighted ? SimpleWeightedGraph() : SimpleGraph()

	add_vertices!(g, 9)

	add_edge!(g, 1, 2)
	add_edge!(g, 2, 3)
	add_edge!(g, 3, 4)
	add_edge!(g, 1, 3)
	add_edge!(g, 3, 5)
	add_edge!(g, 4, 5)
	add_edge!(g, 5, 8)
	add_edge!(g, 5, 7)
	add_edge!(g, 6, 7)
	add_edge!(g, 7, 8)
	add_edge!(g, 7, 9)

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
