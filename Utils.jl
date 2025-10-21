module SDNUtils

include("utils/SdnGraphUtils.jl")
include("utils/SndLib_Parser.jl")

using Random
using Statistics
using StatsBase
using Printf

using GLMakie

using Graphs
using SimpleWeightedGraphs

using LinearAlgebra

using .SNDlibParser
using .SdnGraphUtils

export read_graph, read_sndgraph, safe_probs, to_placement, project_simplex, attack_graph,
    gen_random_vector, game_outcome

const GenericGraph = Union{SimpleGraph, SimpleWeightedGraph}

"""
    read_graph(filename)

    Read the SNDlib graph given FILENAME.
"""
function read_graph(filename)
    SNDlibParser.parse_sndlib_graph(read(filename, String))
end

function read_sndgraph(filename)
    sndg = read_graph(filename)
    idx2node = [n.id for n in sndg.nodes]
	node2idx = Dict(v => i for (i, v) ∈ pairs(idx2node))

	g = SimpleGraph()
	add_vertices!(g, length(node2idx))

	for e ∈ sndg.edges
		α, β = e.source, e.target
		add_edge!(g, node2idx[α], node2idx[β])
	end

	g_pos = [
		Point2f(sndg.nodes[i].lon, sndg.nodes[i].lat)
		for i ∈ 1:length(idx2node)
	]

	g_labels = idx2node

    (
        g = g,
        positions = g_pos,
        labels = g_labels,
        idx2node = idx2node,
        node2idx = node2idx
    )
	# g_sndgraph = SDNGraph(g, g_pos, [], [], g_labels)
end
"""
Project x onto the probability simplex {p ≥ 0, sum(p)=1}.
Handles NaN/Inf and near-zero sums robustly.
"""
function project_simplex(x::AbstractVector{<:Real}; tol=1e-9)
    @assert all(isfinite, x) "non-finite entry in x"
    n = length(x)
    u = sort(collect(x); rev=true)
    css = cumsum(u)
    rho = findlast(i -> u[i] + (1 - css[i]) / i > 0, 1:n)
    theta = (css[rho] - 1) / rho
    p = max.(x .- theta, 0.0)
    s = sum(p)
    if s ≤ tol
        p .= 1/n
    else
        p ./= s
    end
    p
end

"""
The probability `p` returned by the linear program may have values less than `0`. We squeeze these
values between 0 and 1 and normalize them to ensure that we have safe probability that we
can add to `Distribution`.
"""
function safe_probs(p::AbstractVector; tol=1e-12, onzero=:uniform)
    q = copy(p)
    # squash signed zeros and tiny negatives
    q[abs.(q) .< tol] .= 0.0
    if any(q .< 0)
        # allow only "numerical dust" negatives
        if any(q .< -tol)
            throw(DomainError(p, "contains negative mass < -$tol"))
        end
        q[q .< 0] .= 0.0
    end
    s = sum(q)
    if s == 0
        if onzero === :uniform
            q .= 1/length(q)
        else
            throw(DomainError(p, "all mass was zero after clamping"))
        end
    else
        q ./= s
    end
    q
end

"""
Create a placement of ints `[1, 4, 6]` given some one-encoded vector.
"""
function to_placement(placement_bits)
	sort(Int.(findall(BitVector(round.(placement_bits)))))
end


"""
Generate a `V`-sized vector with `K` ones.
"""
function gen_random_vector(V :: Int, K :: Int)
	v = vcat(ones(Int, K), zeros(Int, V - K))
	shuffle!(v)
	sort(to_placement(v))
end

"""
    attack_graph(
    g :: SimpleGraph,
    attack :: AbstractVector{Int}
)

Attack a graph on the nodes in `attack`. This returns the remaining subgraph and
a node mapping which maps the new nodes to the older nodes: `vmap[new_idx] = old_idx`
"""
function attack_graph(
    g :: GenericGraph,
    attack :: AbstractVector{Int}
)
	non_attacked_nodes = setdiff(1:nv(g), attack)
	induced_subgraph(g, non_attacked_nodes)
end

"""
    game_outcome(
    g :: SimpleGraph,
    s :: AbstractVector{Int},
    a :: AbstractVector{Int}
)

Count the number of surviving nodes given the placement `s` and the
attack `a`.
"""
function game_outcome(
    g :: GenericGraph,
    s :: AbstractVector{Int},
    a :: AbstractVector{Int}
)
	# How many nodes survive attack a given placement s
	attacked_g, vmap = attack_graph(g, a)

	surviving_nodes = filter(
		c -> intersect(s, vmap[c]) != [],
		connected_components(attacked_g)
	)

	sum(length.(surviving_nodes))
end

end
