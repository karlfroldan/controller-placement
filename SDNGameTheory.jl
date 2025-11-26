include("Utils.jl")

using JuMP
using CPLEX
using GLMakie
using GraphMakie
using Random
using Distributions
using Statistics
using Combinatorics
using StatsBase
using Printf

using Graphs
using SimpleWeightedGraphs

using LinearAlgebra
import .SDNUtils

DEFAULT_OPTIMIZER = CPLEX.Optimizer
const SEED = 1729

macro time_only(exp)
    quote
        local _t0 = time_ns()
        $(esc(exp))
        local _t1 = time_ns()
        (_t1 - _t0) / 1e9
    end
end

"Returns (value, elapsed seconds)"
macro time_val(expression)
    quote
        local _t0 = time_ns()
        local _val = $(esc(expression))
        local _t1 = time_ns()

        (_val, (_t1 - _t0) / 1e9)
    end
end

"""
    controller_delays(g :: SimpleWeightedGraph, placement :: Vector{Int})

Returns a vector that contains the delays between the controllers.
"""
function controller_delays(g :: SimpleWeightedGraph, placement :: AbstractVector{Int}) Vector{Float64}
    d_table = precalc_delays(g)

    seen = Set{Tuple{Int, Int}}()
    delays = Float64[]

    for (v, w) in Iterators.product(placement, placement)
        if v == w
            continue
        end

        if (w, v) ∉ seen
            push!(delays, d_table[v, w])
            push!(seen, (v, w))
        end
    end

    delays
end

function weight_matrix(graph :: SimpleWeightedGraph)
	n = nv(graph)
	W = fill(Inf, n, n)

	for e ∈ edges(graph)
		u, v = src(e), dst(e)
		w = weight(e)
		W[u, v] = w
		W[v, u] = w
	end 

	@inbounds for i ∈ 1:n
		W[i, i] = 0.0
	end
	W
end

function precalc_delays(graph :: SimpleWeightedGraph)
	n = nv(graph)
	W = weight_matrix(graph)
	# Delay is represented by the matrix of size VxV where the path from v to w
	# has delay D_{v, w}
	D = fill(Inf, n, n)

	for v ∈ vertices(graph)
		state = dijkstra_shortest_paths(graph, v, W)
		D[v, :] = state.dists
	end
	D
end

function minmax_sc_delay(g :: SimpleWeightedGraph, M :: Int, bsc = 0.0, bcc = 0.0)
    m = Model(DEFAULT_OPTIMIZER)
    set_silent(m)

    delay_table = precalc_delays(g)
    V = nv(g)
    # @show delay_table

    @variable(m, D)
    # Controller placements
    @variable(m, s[1:V], Bin)
    # z_vw = 1 iff switch in v is assigned to controller w
    # @variable(m, z[
    Ws = Dict(v => sort([w for w in vertices(g) if delay_table[v, w] ≤ bsc]) for v in vertices(g))

    z = [@variable(m, [Ws[v]], Bin) for v in vertices(g)]

    # (A1b)
    @constraint(m, sum(s) == M) # M controllers required

    # (A1c)
    U = Vector{Tuple{Int, Int}}()
    for (v, w) in Iterators.product(vertices(g), vertices(g))
        if delay_table[v, w] > bcc && (w, v) ∉ U
            push!(U, (v, w))
            @constraint(m, s[v] + s[w] ≤ 1)
        end
    end

    # (A1d)
    for (k, w) in Ws
        @constraint(m, sum(z[k]) == 1)
    end

    for v in vertices(g)
        ws = setdiff(Ws[v], [v])
        
        # (A1e)
        @constraint(m, z[v][ws] .≤ s[ws])
        # (A1f)
        @constraint(m, z[v][v] == s[v])

        # (A1g)
        ds = delay_table[v, Ws[v]]
        @constraint(m, D ≥ sum(ds .* z[v]))
    end

    @objective(m, Min, D)

    optimize!(m)
    # @show objective_value(m)
    
    # m, D
    value(D)
end

function cpop(
	g :: AbstractGraph, 
	num_controllers :: Int,
	attacks :: Vector{Vector{Int}};
	solver = DEFAULT_OPTIMIZER,
)
	m = Model(solver)
	V = nv(g)
	attack_size = length(attacks)

	set_silent(m)

	# Controller placement variables
	@variable(m, s[1:V], Bin)
	# Variable to determine whether the node survives attack
	@variable(m, y[1:V, 1:attack_size], Bin)
	# Value to maximize
	@variable(m, Y ≥ 0)

	@objective(m, Max, Y)

	# Constraint on the number of controllers
	@constraint(m, sum(s) == num_controllers)

	# All nodes that were attacked should not have survived that attack
	for (i, a) ∈ enumerate(attacks)
		@constraint(m, [v in a], y[v, i] == 0)
	end

	# If the component c does not contain any controller, then every y_{v, a} = 0
	# First, let's get all the components resulting from a given attack a
	for (i, a) ∈ enumerate(attacks)
		g_attacked, vmap = SDNUtils.attack_graph(g, a)
		C = connected_components(g_attacked)

		@constraint(m, [c in C], sum(y[vmap[c], i]) ≤ length(c) * sum(s[vmap[c]]))
	end

	# Last constraint : For every attack, maximize the number of surviving nodes
	# Above Y, in the best case
	for i ∈ eachindex(attacks)
		@constraint(m, Y ≤ sum(y[:, i]))
	end

	m, s
end

function cpop_with_delays(
	g :: SimpleWeightedGraph, 
	num_controllers :: Int,
	attacks :: Vector{Vector{Int}},
    bsc :: Float64,
    bcc :: Float64;
	solver = DEFAULT_OPTIMIZER,
)
	m = Model(solver)
	V = nv(g)
	attack_size = length(attacks)

	set_silent(m)

    delay_table = precalc_delays(g)

	# Controller placement variables
	@variable(m, s[1:V], Bin)
	# Variable to determine whether the node survives attack
	@variable(m, y[1:V, 1:attack_size], Bin)
	# Value to maximize
	@variable(m, Y ≥ 0)

	@objective(m, Max, Y)

	# Constraint on the number of controllers
	@constraint(m, sum(s) == num_controllers)

	# All nodes that were attacked should not have survived that attack
	for (i, a) ∈ enumerate(attacks)
		@constraint(m, [v in a], y[v, i] == 0)
	end

    # Set of nodes that for every v in vertex set,
    # satisfies the BSC constraint
    for v in vertices(g)
        d = delay_table[v, :] .< bsc
        @constraint(m, sum(d .* s) ≥ 1)
    end

    # Get all node pairs that violate BCC
    U = Vector{Tuple{Int, Int}}()
    for (v, w) in Iterators.product(vertices(g), vertices(g))
        if delay_table[v, w] > bcc && (w, v) ∉ U
            push!(U, (v, w))
            @constraint(m, s[v] + s[w] ≤ 1)
        end
    end

	# If the component c does not contain any controller, then every y_{v, a} = 0
	# First, let's get all the components resulting from a given attack a
	for (i, a) ∈ enumerate(attacks)
		g_attacked, vmap = SDNUtils.attack_graph(g, a)
		C = connected_components(g_attacked)

		@constraint(m, [c in C], sum(y[vmap[c], i]) ≤ length(c) * sum(s[vmap[c]]))
	end

	# Last constraint : For every attack, maximize the number of surviving nodes
	# Above Y, in the best case
	for i ∈ eachindex(attacks)
		@constraint(m, Y ≤ sum(y[:, i]))
	end

	m, s
end

function naop(
	g :: AbstractGraph,
	K :: Int,
	controller_placements :: Vector{Vector{Int}};
	solver = DEFAULT_OPTIMIZER,
)
	m = Model(solver)

	set_silent(m)

	V = nv(g)
	E = ne(g)

	# How many nodes to attack
	@variable(m, a[1:V], Bin)
	# Is link down?
	@variable(m, t[1:E], Bin)
	# 1 if node v is up after an attack on placement s
	@variable(m, z[1:length(controller_placements), 1:V], Bin)
	
	# Objective variable
	@variable(m, Z >= 0)
	
	@objective(m, Min, Z)

	# We need to attack K nodes in total.
	@constraint(m, sum(a) == K)

	# Constraint that ensures that when a node is attacked,
	# all links to it are down
	for v ∈ 1:V
		for (i, e) ∈ enumerate(edges(g))
			α, β = src(e), dst(e)
			if α == v || β == v
				# Incident
				@constraint(m, t[i] ≥ a[v])
			end
		end
	end

	# If link is up t_e = 0, then both sides of the link should be up
	# Otherwise, RHS can either be 1 or 2
	for (i, e) ∈ enumerate(edges(g))
		α, β = src(e), dst(e)
		@constraint(m, t[i] ≤ a[α] + a[β])
	end

	
	for (i, s) ∈ enumerate(controller_placements)
		# If a node v is attacked, then the link should be down for every placement
		@constraint(m, z[i, :] .≤ 1 .- a)
		# If a controller node v is not attacked, then it should survive 
		# and equal to 1
		@constraint(m, z[i, s] .≥ 1 .- a[s])
	end

	# Given a placement s, if the link is up t_e == 0, then both nodes
	# at either end are simultaneously up and vice versa.
	for (i, e) ∈ enumerate(edges(g))
		α, β = src(e), dst(e)
		@constraint(m, [j in eachindex(controller_placements)], 
					z[j, α] ≥ z[j, β] - t[i])
		@constraint(m, [j in eachindex(controller_placements)], 
					z[j, β] ≥ z[j, α] - t[i])
	end

	# The number of surviving nodes should be lower than Z
	@constraint(m, [i in eachindex(controller_placements)],
			   Z ≥ sum(z[i, :]))
	
	m, a
end

function controller_placement_optimization(
    g :: AbstractGraph, M :: Int, K :: Int;
    TOL :: Float64 = 1e-6, has_delay = false, bcc = 0.0, bsc = 0.0
    )
    if has_delay
        @assert(bcc != 0.0)
        @assert(bsc != 0.0)
    end

	attack_set = Vector{Int}[]
	Y_star = Float64(nv(g))
	Z_star = Float64(nv(g))

	# STEP 0: Generate a random M-node controller placement s*
	s_star = SDNUtils.gen_random_vector(nv(g), M)
	n_iterations = 0
    naop_times = []
    cpop_times = []

	# STEP 1: Solve A[K, {s*}] to get the worst attack a*. If Z* >= Y* go to step 3
	while true
		# Solve NAOP for the current best placement
		placements = [s_star]
		naop_model, a_star = naop(g, K, placements)

		elapsed = @time_only(optimize!(naop_model))

		assert_is_solved_and_feasible(naop_model)

        push!(naop_times, elapsed)
		
		a_star = SDNUtils.to_placement(value.(a_star))
		Z_star = objective_value(naop_model)

		if Z_star >= Y_star - TOL
			break
		end

		# STEP 2: A = A union {a*}. Then solve P[M, A] to get placement s*. Repeat
		if a_star ∉ attack_set
			push!(attack_set, a_star)
		else
            error("Attack [$a_star] is already in the attack set")
		end

        if has_delay
            cpop_model, s_star = cpop_with_delays(g, M, attack_set, bsc, bcc)
        else
		    cpop_model, s_star = cpop(g, M, attack_set)
        end
        
		elapsed = @time_only(optimize!(cpop_model))

		assert_is_solved_and_feasible(cpop_model)

		s_star = SDNUtils.to_placement(value.(s_star))
		Y_star = objective_value(cpop_model)
		n_iterations += 1

        push!(cpop_times, elapsed)
	end

	s_star, Y_star
end

function attack_optimization(
    g :: AbstractGraph, M :: Int, K :: Int;
    TOL :: Float64 = 1e-6, has_delay = true, bsc = 0.0, bcc = 0.0
    )
    if has_delay
        @assert(bsc != 0.0)
        @assert(bcc != 0.0)
    end

	placement_set = Vector{Int}[]
	Z_star = 0.0
	Y_star = 0.0

	# STEP 0: Generate a random K-node attack a*
	a_star = vcat(ones(Int, K), zeros(Int, nv(g) - K))
	shuffle!(a_star)
	a_star = SDNUtils.to_placement(a_star)
	n_iterations = 0

	cpop_times = []
	naop_times = []

	# STEP 1: Solve P[M, {a*}] to get the best placement s*. If Y* <= Z*, STEP 4
	while true
		# Solve CPOP for the current best placement
		attack_set = [a_star]
        if has_delay
            cpop_model, s_star = cpop_with_delays(g, M, attack_set, bsc, bcc)
        else
		    cpop_model, s_star = cpop(g, M, attack_set)
        end

		elapsed = @time_only(optimize!(cpop_model))
		push!(cpop_times, elapsed)

		is_solved_and_feasible(cpop_model)

		s_star = SDNUtils.to_placement(value.(s_star))
		Y_star = objective_value(cpop_model)

		if Y_star <= Z_star + TOL
			break
		end

		# STEP 2: S = S union {s*}. Then solve A[K, S] to get attack a*. Repeat
		if s_star ∉ placement_set
			push!(placement_set, s_star)
		else
            error("placement [$s_star] is already in the placement set.")
		end

		naop_model, a_star = naop(g, K, placement_set)
		elapsed = @time_only(optimize!(naop_model))
		push!(naop_times, elapsed)

		is_solved_and_feasible(naop_model)

		a_star = SDNUtils.to_placement(value.(a_star))
		Z_star = objective_value(naop_model)

		n_iterations += 1
	end

	# STEP 3: STOP. Current Placement s* is an optimal solution
	
	a_star, Z_star
end


function master_placement(
    g :: AbstractGraph, placements :: Vector{Vector{Int}}, attacks :: Vector{Vector{Int}}
)
	# P[S, A] = Find optimal decision probability decision for the defender over the
	# set of all placements and attacks

	m = Model(DEFAULT_OPTIMIZER)
	set_silent(m)

	n_placements = length(placements)

	@variable(m, y)
	@variable(m, q[1:n_placements] ≥ 0)

	@objective(m, Max, y)

    @constraint(m, sum(q) == 1)

	# Game outcome constraints. Their dual is p
	p_dual = ConstraintRef[]
	for a ∈ attacks
		vs = [SDNUtils.game_outcome(g, s, a) for s ∈ placements]
		c = @constraint(m, y ≤ sum(vs .* q))
		push!(p_dual, c)
	end
	
	m, q, p_dual
end

function pricing_placement_aux(
	g :: AbstractGraph,
	M :: Int,
	attack_set :: Vector{Vector{Int}},
	p_star :: Vector{Float64}
)
	# M - M controllers required.
	# p_star - Current attacker's probability distribution
	
	V = nv(g)
	m = Model(DEFAULT_OPTIMIZER)
	set_silent(m)
	attack_len = length(attack_set)

	# Total number of nodes surviving attack a
	@variable(m, Y[1:attack_len] ≥ 0, Int)
	# Node is chosen to host a controller
	@variable(m, s[1:V], Bin)
	# y[v, a] = 1 <=> node v survives attack a when constructed placement s is used.
	@variable(m, y[1:V, 1:attack_len], Bin)

	@objective(m, Max, sum(p_star .* Y))

	# (24b) There should be equal to M controllers
	@constraint(m, sum(s) == M)

	for (i, a) ∈ enumerate(attack_set)
		# (24c) All attacked nodes are 0
		@constraint(m, [v in a], y[v, i] == 0)

		# (24d) Try to ensure that each component has a surviving controller
		attacked_g, vmap = SDNUtils.attack_graph(g, a)
		for c ∈ connected_components(attacked_g)
			c′ = vmap[c]
			@constraint(m, sum(y[c′, i]) ≤ length(c′) * sum(s[c′]))
		end

		# (24e) Make y and Y correspond to each other 
		@constraint(m, Y[i] == sum(y[:, i]))
	end

	(
        model = m,
        s = s,
        y = y,
        Y = y
    )
end


function pricing_placement(
    g :: SimpleWeightedGraph,
    M :: Int,
    attack_set :: Vector{Vector{Int}},
    p_star :: Vector{Float64},
    delay_table :: Matrix{Float64},
    bcc :: Float64,
    bsc :: Float64
    )
    # V = nv(g)
    
    out = pricing_placement_aux(g, M, attack_set, p_star)
    m = out.model
    s = out.s
    # Y = out.Y
    # y = out.y

    # Set of nodes that for every v in vertex set,
    # satisfies the BSC constraint
    for v in vertices(g)
        d = delay_table[v, :] .< bsc
        @constraint(m, sum(d .* s) ≥ 1)
    end

    # Get all node pairs that violate BCC
    U = Vector{Tuple{Int, Int}}()
    for (v, w) in Iterators.product(vertices(g), vertices(g))
        if delay_table[v, w] > bcc && (w, v) ∉ U
            push!(U, (v, w))
            @constraint(m, s[v] + s[w] ≤ 1)
        end
    end

    m, s
end

function pricing_placement(
	g :: AbstractGraph,
	M :: Int,
	attack_set :: Vector{Vector{Int}},
	p_star :: Vector{Float64}
)
    out = pricing_placement_aux(g, M, attack_set, p_star)
    out.model, out.s
end

function pricing_attack(
	g :: AbstractGraph,
	K :: Int,
	controller_set :: Vector{Vector{Int}},
	q_star :: Vector{Float64}
)
	m = Model(DEFAULT_OPTIMIZER)
	set_silent(m)

	# Consider only the controller placements with q > 0
	S = [c for (q, c) in zip(q_star, controller_set) if q > 0]
	q_star_prime = [q for q in q_star if q > 0]

	V = nv(g)

	# How many nodes survive given placement s
	@variable(m, F[1:length(S)] ≥ 0, Int)
	# Attack placement variable
	@variable(m, a[1:V], Bin)
	# Does node v survive the attack given placement s?
	@variable(m, z[1:V, 1:length(S)], Bin)

	@objective(m, Min, sum(q_star_prime .* F))

	# (25b) Attack exactly K nodes
	@constraint(m, sum(a) == K)

	# (25d) If controller β(e) does not survive the attack, then
	#       its neighbors should also not survive the attacks.
	# (25e) Same but reversed.
	for e ∈ edges(g)
		α, β = src(e), dst(e)
		@constraint(m, z[β, :] .≥ z[α, :] .- a[β])
		@constraint(m, z[α, :] .≥ z[β, :] .- a[α])
	end

	
	for (i, s) ∈ enumerate(S)
		# (25c) Ensure that if node v is attacked, then it shouldn't survive
		@constraint(m, z[s, i] .== 1 .- a[s])

		# (25f) Synchronize z and F
		@constraint(m, F[i] == sum(z[:, i]))
	end
	
	m, a
end

function find_feasible_controller_placement(g :: AbstractGraph, M :: Int)
    SDNUtils.gen_random_vector(nv(g), M)
end

"""
    switch2controller_delay(delay_table :: Matrix{Float64}, s :: Vector{Float64})

Calculate the switch-to-controller given switch `switch`
"""
function get_switch2controller_delay(delay_table :: Matrix{Float64}, s :: Vector{Int}, switch :: Int)
    minimum(delay_table[switch, s])
end

function get_controller2controler_delay(delay_table :: Matrix{Float64}, s :: Vector{Int}, controller :: Int)
    @assert(controller in s)

    minimum(delay_table[controller, setdiff(s, controller)])
end

function switch2controller_delays(delay_table :: Matrix{Float64}, s :: Vector{Int})
    [get_switch2controller_delay(delay_table, s, i) for i in 1:size(delay_table, 1)]
end

function controller2controller_delays(delay_table :: Matrix{Float64}, s :: Vector{Int})
    [get_controller2controler_delay(delay_table, s, i) for i in s]
end

function find_feasible_controller_placement(g :: SimpleWeightedGraph, M :: Int, bcc :: Float64, bsc :: Float64)
    m = Model(DEFAULT_OPTIMIZER)
    set_silent(m)

    delay_table = precalc_delays(g)

    V = nv(g)

    @variable(m, s[1:V], Bin)
    @constraint(m, sum(s) == M)
    
    # Set of nodes that for every v in vertex set,
    # satisfies the BSC constraint
    for v in vertices(g)
        d = delay_table[v, :] .< bsc
        @constraint(m, sum(d .* s) ≥ 1)
    end

    # Get all node pairs that violate BCC
    U = Vector{Tuple{Int, Int}}()
    for (v, w) in Iterators.product(vertices(g), vertices(g))
        if delay_table[v, w] > bcc && (w, v) ∉ U
            push!(U, (v, w))
            @constraint(m, s[v] + s[w] ≤ 1)
        end
    end
    w = rand(V)                 # in (0,1)
    @objective(m, Min, sum(w[i]*s[i] for i in 1:V))
    # @objective(m, Min, 0)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        return SDNUtils.to_placement(value.(s))
    else
        return nothing
    end
end

function mixed_strategy_algorithm(
	g :: AbstractGraph, 
	M :: Int, K :: Int;
    TOL = 1e-8,
    bcc :: Float64 = 0.0,
    bsc :: Float64 = 0.0,
    has_delays = false
)
    if has_delays
        @assert(bsc != 0.0)
        @assert(bcc != 0.0)
        delay_table = precalc_delays(g)
    end
	
	# Random attack and random placement of size M and K
    if has_delays
        s = find_feasible_controller_placement(g, M, bcc, bsc)
    else
        s = find_feasible_controller_placement(g, M)
    end
    
	a = SDNUtils.gen_random_vector(nv(g), K)

    @assert(!isnothing(s))

	placements = [s]
	attacks = [a]

	master, q, shadow_prices = master_placement(g, placements, attacks)
	optimize!(master)
	assert_is_solved_and_feasible(master)

	# Obviously, both will be [1.0] at first
	q = value.(q)
	p = SDNUtils.project_simplex(-dual.(shadow_prices))

	x_star = dual_objective_value(master)
	y_star = objective_value(master)
    
    n_iterations = 0

    has_updated = true

    time_placement = Float64[]
    time_attack = Float64[]
    time_master = Float64[]

    stat_v_star = Float64[]
    stat_expected_value = Float64[]

	while has_updated
        has_updated = false
		
		# Solve the placement generation problem to get a new placement s'
        if has_delays
            cp, s = pricing_placement(g, M, attacks, p, delay_table, bcc, bsc)
        else
		    cp, s = pricing_placement(g, M, attacks, p)
        end

		elapsed = @time_only(optimize!(cp))
        push!(time_placement, elapsed)
        
		assert_is_solved_and_feasible(cp)
		s = SDNUtils.to_placement(value.(s))
		# Generated a new placement s
		outcomes = [SDNUtils.game_outcome(g, s, a′) for a′ ∈ attacks]
		expected = outcomes ⋅ p

		if expected > x_star - TOL
			if s ∉ placements
				has_updated = true
				push!(placements, s)
			end
		end

		# Then solve the master problem again.
		master, q, shadow_prices = master_placement(g, placements, attacks)
		time_master_elapsed = @time_only(optimize!(master))
        
		assert_is_solved_and_feasible(master)

		q = value.(q)
		p = SDNUtils.project_simplex(-dual.(shadow_prices))

		x_star = dual_objective_value(master)
		y_star = objective_value(master)

        @assert(x_star >= y_star - TOL)

		# Then generate a new attack
		na, a = pricing_attack(g, K, placements, q)
		elapsed = @time_only(optimize!(na))
        push!(time_attack, elapsed)

		assert_is_solved_and_feasible(na)

		a = SDNUtils.to_placement(value.(a))
		outcomes = [SDNUtils.game_outcome(g, s′, a) for s′ ∈ placements]
		expected = outcomes ⋅ q
		
		if y_star > expected - TOL
			if a ∉ attacks
				has_updated = true
				push!(attacks, a)
			end
		end

		master, q, shadow_prices = master_placement(g, placements, attacks)
		time_master_elapsed += @time_only(optimize!(master))
		assert_is_solved_and_feasible(master)
        
        push!(time_master, time_master_elapsed)

		q = value.(q)
		p = SDNUtils.project_simplex(-dual.(shadow_prices))
		x_star = dual_objective_value(master)
		y_star = objective_value(master)

        @assert(x_star >= y_star - TOL)

        push!(stat_v_star, x_star)
        push!(stat_expected_value, expected)

        n_iterations += 1
	end

    (
        p_star = p,
        q_star = q,
        V_star = x_star,
        attacks = attacks,
        placements = placements,
        stats_v_star = stat_v_star,
        stats_expected_values = stat_expected_value,
        stats_time_placement = time_placement,
        stats_time_attack = time_attack,
        stats_time_master = time_master
    )
end

function outcome_metrics_mixed(g, M::Int, K::Int)
    V = nv(g)
    0 ≤ M ≤ V || throw(ArgumentError("M must be in 0..$V"))
    0 ≤ K ≤ V || throw(ArgumentError("K must be in 0..$V"))

    ms = mixed_strategy_algorithm(g, M, K)
    V_star = ms.V_star
	return V_star
end

function outcome_metrics_pure(g :: AbstractGraph, M::Int, K::Int)
    V = nv(g)
    @assert 0 ≤ M ≤ V && 0 ≤ K ≤ V
	_, v_maxmin = controller_placement_optimization(g, M, K)
	_, v_minmax = attack_optimization(g, M, K)

    return v_minmax, v_maxmin
end

function outcome_table(g :: AbstractGraph, M_max :: Int, K_max :: Int; K_min :: Int = 2, full = false)
	v_minmax_table = Matrix{Int}(undef, M_max, K_max)
	v_maxmin_table = Matrix{Int}(undef, M_max, K_max)
	v_star_table = Matrix{Float64}(undef, M_max, K_max)
	
	for m ∈ 1:M_max
		for k ∈ K_min:K_max
            @printf("---- Performing K = %d, M = %d ---- \n", k, m)
			v_star = outcome_metrics_mixed(g, m, k)
			v_star_table[m, k] = v_star
			(v_minmax, v_maxmin) = outcome_metrics_pure(g, m, k)
			v_minmax_table[m, k] = v_minmax
			v_maxmin_table[m, k] = v_maxmin
            @printf("\n\n")
		end
	end

    if full
        full_matrix = Matrix{Tuple{Int, Float64, Int}}(undef, M_max, K_max)
        for m in 1:M_max
            for k in K_min:K_max
                full_matrix[m, k] = (v_maxmin_table[m, k], v_star_table[m, k], v_minmax_table[m, k])
            end
        end

        return full_matrix[:, K_min:K_max]
    else
	    return (
            minmax = v_minmax_table[:, K_min:K_max],
	        maxmin = v_maxmin_table[:, K_min:K_max],
	        v_star = v_star_table
        )
    end
end

"""
    sample_strategy(action :: AbstractVector{Vector{Int64}}, distribution :: AbstractVector{Float64})

Given a list of strategies and a corresponding probability distribution, sample a strategy.
"""
function sample_strategy(action :: AbstractVector{Vector{Int64}}, distribution :: AbstractVector{Float64})
    TOL = 1e-9
    @assert(sum(distribution) >= 1.0 - TOL)
    @assert(length(action) == length(distribution))

    dist = Distributions.Categorical(distribution)
    action[rand(dist)]
end

function plot_vstar_and_expected(v_star :: AbstractVector{Float64}, expected :: AbstractVector{Float64})
    @assert(length(v_star) == length(expected))

    n = length(v_star)

    f = Figure()
    ax = Axis(f[1, 1])
    lines!(ax, 1:n, v_star)
    lines!(ax, 1:n, expected)

    f
end

function general_stats(g :: AbstractGraph, s :: Vector{Int})
    delays = precalc_delays(g)
    sc = switch2controller_delays(delays, s)
    cc = controller_delays(g, s)

    (
        mean_sc = mean(sc),
        max_sc = maximum(sc),
        mean_cc = mean(cc),
        max_cc = maximum(cc)
    )
end

const ExpectedValueNothing = Union{Float64, Symbol}

function safe_mixed_strategy(args...; kwargs...)
    try
        mixed_strategy_algorithm(args...; kwargs...)
    catch
        :infeasible
    end
end

# For simple.g
# mixed_strategies_stats(simple.g, 270, 350, 450, 250, 3, 3)
# It's good practice to define the Union type for clarity
# const ExpectedValueNothing = Union{Float64, Symbol}

struct MixedStrategyStatistic
    payoff :: Union{Float64, Symbol}
    # The number of non-zero placements
    n_placements :: Int
    # The number of non-zero attacks
    n_attacks :: Int
    # Runtime in seconds
    exec_time :: Float64
    placement_entropy :: Float64
    attack_entropy :: Float64
    iterations :: Int
end

# Helper for the minmax calculation, mirroring safe_mixed_strategy
function safe_minmax_sc_delay(g, m, bsc, bcc; TOL=1e-9)
    try
        minmax_sc_delay(g, m, bsc, bcc) + TOL
    catch
        :infeasible
    end
end

infeasible_mixed_strategy_statistic(time) = MixedStrategyStatistic(:infeasible, 0, 0, time, 0, 0, 0)

function entropy(arr :: Vector{Float64}; base=2)
    # @assert(arr .>= 0.0 && arr .<= 1.0)

    filtered = filter(x -> x > 1e-9, arr)

    if length(filtered) == 0
        return 0.0
    end

    - sum(filtered .|> p -> p * log(base, p))
end

function mixed_strategies_stats(
    g::SimpleWeightedGraph,
    bcc_low::Float64,
    bcc_high::Float64,
    bsc::Float64,
    m_high::Int,
    k_high::Int;
    k_low::Int = 1,
    m_low::Int = 1,
    TOL::Float64 = 1e-9,
    )
    
    n_non_zero(arr) = length(filter(x -> x > 1e-9, arr))
    n_iterations(stats) = length(stats.stats_v_star)
    safe_probs = SDNUtils.safe_probs

    new_statistic(stat, time) =
        MixedStrategyStatistic(
            stat.V_star,
            n_non_zero(safe_probs(stat.q_star)),
            n_non_zero(safe_probs(stat.p_star)),
            time,
            entropy(safe_probs(stat.q_star)),
            entropy(safe_probs(stat.p_star)),
            n_iterations(stat),
        )

    return [
        begin
            local low_stats, high_stats, no_delay_stats
            @show k, m
            # 1. Calculate bsc values safely
            bsc_low = safe_minmax_sc_delay(g, m, bsc, bcc_low; TOL=TOL)
            bsc_high = safe_minmax_sc_delay(g, m, bsc, bcc_high; TOL=TOL)

            # 2. Check for infeasibility before calculating the next steps
            low, low_time = @time_val((bsc_low == :infeasible) ? :infeasible :
                safe_mixed_strategy(g, m, k; has_delays=true, bcc=bcc_low, bsc=bsc_low))

            if low == :infeasible
                low_stats = infeasible_mixed_strategy_statistic(low_time)
            else
                low_stats = new_statistic(low, low_time)
            end

            high, high_time = @time_val((bsc_high == :infeasible) ? :infeasible :
                safe_mixed_strategy(g, m, k; has_delays=true, bcc=bcc_high, bsc=bsc_high))

            if high == :infeasible
                high_stats = infeasible_mixed_strategy_statistic(high_time)
            else
                high_stats = new_statistic(high, high_time)
            end

            # This one is independent and can be calculated directly
            no_delay, no_delay_time = @time_val(safe_mixed_strategy(g, m, k; has_delays=false))
            if no_delay == :infeasible
                no_delay_stats = infeasible_mixed_strategy_statistic(no_delay_time)
            else
                no_delay_stats = new_statistic(no_delay, no_delay_time)
            end

            @printf("[Exec Times] low = %.2f, no_delay = %.2f, high = %.2f\n",
                    low_time, no_delay_time, high_time)

            # 3. Assemble the final tuple
            (low_stats, no_delay_stats, high_stats)
            # (bsc_low, low, no_delay, high, bsc_high)
        end
        # The comprehension builds the correctly sized matrix directly
        for m in m_low:m_high, k in k_low:k_high
    ]
end




"""
_format_latex_value(stats::MixedStrategyStatistic, column::Symbol, decimal_places::Int)

Internal helper function to format the values from the MixedStrategyStatistic struct.
Handles :infeasible symbols and formats floats to a fixed decimal place.
"""
function _format_latex_value(stats::MixedStrategyStatistic, column::Symbol, decimal_places::Int)
    val = getfield(stats, column)
    
    if val == :infeasible
        return raw"\textemdash"
    elseif val isa AbstractFloat
        return @sprintf "%.*f" decimal_places val
    elseif val isa Int
        return string(val)
    else
        return string(val)
    end
end

"""
    stats_to_latex(
        mat :: Matrix{Tuple{MixedStrategyStatistic, MixedStrategyStatistic, MixedStrategyStatistic}}, 
        filename :: AbstractString, 
        table_title :: AbstractString, 
        column :: Symbol; 
        k_start=1, 
        m_start=1, 
        decimal_places=2,
        k_cols_per_table=3,
        row_labels = ("(BCC=1500)", "No Delay", "(BCC=2000)")
    )

Generates a complete, self-contained LaTeX table environment from the mixed strategy 
statistics matrix and writes it to `filename`.

If the number of K columns exceeds `k_cols_per_table`, it will generate
multiple `tabular` environments (chunks) within the same `table` environment.

The `column` argument is a `Symbol` specifying which field from the 
`MixedStrategyStatistic` struct to display (e.g., `:payoff`, `:exec_time`).
"""
function stats_to_latex(
    mat :: Matrix{Tuple{MixedStrategyStatistic, MixedStrategyStatistic, MixedStrategyStatistic}}, 
    filename :: AbstractString, 
    table_title :: AbstractString,
    column :: Symbol; 
    k_start=1, 
    m_start=1, 
    decimal_places=2,
    k_cols_per_table=3,
    row_labels = ("Tight Delay (BCC=1500)", "No Delay", "Loose Delay (BCC=2000)")
)
    
    open(filename, "w") do file
        n_rows, n_cols = size(mat) # n_rows = M count, n_cols = K count
        
        # Determine how many table chunks to create
        num_blocks = ceil(Int, n_cols / k_cols_per_table)
        
        # --- Table Preamble ---
        write(file, raw"% Add \usepackage{booktabs}, \usepackage{multirow}, \usepackage{amsmath} to your LaTeX preamble" * "\n")
        write(file, raw"\begin{table}[!htbp]" * "\n")
        write(file, raw"\centering" * "\n")
        
        write(file, "\\caption{$(table_title)}\n")
        label_name = "tab:" * lowercase(replace(table_title, r"[^a-zA-Z0-9\s]" => "", " " => "-"))
        write(file, "\\label{$(label_name)}\n")
        
        write(file, raw"\small" * "\n\n")

        for block_idx in 1:num_blocks
            k_start_idx = (block_idx - 1) * k_cols_per_table + 1
            k_end_idx = min(block_idx * k_cols_per_table, n_cols)
            current_k_indices = k_start_idx:k_end_idx
            num_current_cols = length(current_k_indices)

            col_specifiers = repeat("c", num_current_cols)
            col_format = "@{} l c $(col_specifiers) @{}"
            write(file, "\\begin{tabular}{$(col_format)}\n")
            write(file, raw"\toprule" * "\n")

            write(file, "\\textbf{Delay Setting} & \\textbf{M} & \\multicolumn{$(num_current_cols)}{c}{\\textbf{K (Attacks)}} \\\\\n")
            
            write(file, " & & ") # Skip first two columns
            k_headers = ["\\textbf{$(k_start + k_idx - 1)}" for k_idx in current_k_indices]
            write(file, join(k_headers, " & ") * " \\\\\n")
            
            write(file, "\\cmidrule(l){3-$(2 + num_current_cols)}\n")

            for m_idx in 1:n_rows
                current_m = m_start + m_idx - 1
                
                # --- Sub-row 1: (low_stats) ---
                write(file, "$(row_labels[1]) & \\multirow{3}{*}{$(current_m)} & ")
                data_row_1 = [
                    _format_latex_value(mat[m_idx, k_idx][1], column, decimal_places) 
                    for k_idx in current_k_indices
                ]
                write(file, join(data_row_1, " & ") * " \\\\\n")

                # --- Sub-row 2: (no_delay_stats) ---
                write(file, "$(row_labels[2]) & & ") # M-column is blank
                data_row_2 = [
                    _format_latex_value(mat[m_idx, k_idx][2], column, decimal_places) 
                    for k_idx in current_k_indices
                ]
                write(file, join(data_row_2, " & ") * " \\\\\n")

                # --- Sub-row 3: (high_stats) ---
                write(file, "$(row_labels[3]) & & ") # M-column is blank
                data_row_3 = [
                    _format_latex_value(mat[m_idx, k_idx][3], column, decimal_places) 
                    for k_idx in current_k_indices
                ]
                write(file, join(data_row_3, " & ") * " \\\\\n")

                # --- Separator ---
                if m_idx < n_rows
                    write(file, raw"\midrule" * "\n")
                end
            end

            # --- Table Footer for this chunk ---
            write(file, raw"\bottomrule" * "\n")
            write(file, raw"\end{tabular}" * "\n")

            # Add space if there is another chunk
            if block_idx < num_blocks
                write(file, raw"\vspace{1em}" * "\n\n") # 1em vertical space
            end
        end

        # --- Close the main table environment ---
        write(file, raw"\end{table}" * "\n")
    end # close file
    
    println("Successfully generated full LaTeX table at $(filename)")
end


"""
    stats_to_latex_dual(
        mat :: Matrix{Tuple{MixedStrategyStatistic, MixedStrategyStatistic, MixedStrategyStatistic}}, 
        filename :: AbstractString, 
        table_title :: AbstractString, 
        column1 :: Symbol, 
        col1_label :: String,
        column2 :: Symbol,
        col2_label :: String;
        k_start=1, 
        m_start=1, 
        decimal_places=2,
        k_cols_per_table=3,
        row_labels = ("(BCC=1500)", "No Delay", "(BCC=2000)")
    )

Generates a complete, self-contained LaTeX table environment for dual (side-by-side)
metrics and writes it to `filename`.

If the number of K columns exceeds `k_cols_per_table`, it will generate
multiple `tabular` environments (chunks) within the same `table` environment.
"""
function stats_to_latex_dual(
    mat :: Matrix{Tuple{MixedStrategyStatistic, MixedStrategyStatistic, MixedStrategyStatistic}}, 
    filename :: AbstractString, 
    table_title :: AbstractString,
    column1 :: Symbol, 
    col1_label :: String,
    column2 :: Symbol,
    col2_label :: String;
    k_start=1, 
    m_start=1, 
    decimal_places=2,
    k_cols_per_table=3,
    row_labels = ("Tight Delay (BCC=1500)", "No Delay", "Loose Delay (BCC=2000)")
)
    
    open(filename, "w") do file
        n_rows, n_cols = size(mat) # n_rows = M count, n_cols = K count

        num_blocks = ceil(Int, n_cols / k_cols_per_table)
        
        write(file, raw"% Add \usepackage{booktabs}, \usepackage{multirow}, \usepackage{amsmath} to your LaTeX preamble" * "\n")
        write(file, raw"\begin{table}[!htbp]" * "\n")
        write(file, raw"\centering" * "\n")
        
        write(file, "\\caption{$(table_title)}\n")
        label_name = "tab:" * lowercase(replace(table_title, r"[^a-zA-Z0-9\s]" => "", " " => "-"))
        write(file, "\\label{$(label_name)}\n")
        
        write(file, raw"\small" * "\n\n")

        # --- Main loop to create table chunks ---
        for block_idx in 1:num_blocks
            # Calculate the K column indices for this chunk
            k_start_idx = (block_idx - 1) * k_cols_per_table + 1
            k_end_idx = min(block_idx * k_cols_per_table, n_cols)
            current_k_indices = k_start_idx:k_end_idx
            num_current_cols = length(current_k_indices)

            col_specifiers = repeat("c", num_current_cols * 2) 
            col_format = "@{} l c $(col_specifiers) @{}"
            write(file, "\\begin{tabular}{$(col_format)}\n")
            write(file, raw"\toprule" * "\n")

            # --- Column Headers ---
            # Row 1: Delay Setting | M | K = 2 | K = 3 | K = 4 ...
            write(file, "\\textbf{Delay Setting} & \\textbf{M} & ")
            k_headers = [
                "\\multicolumn{2}{c}{\\textbf{K = $(k_start + k_idx - 1)}}" 
                for k_idx in current_k_indices
            ]
            write(file, join(k_headers, " & ") * " \\\\\n")
            
            write(file, " & & ") # Skip first two columns
            sub_headers = repeat(["\\textbf{$(col1_label)}", "\\textbf{$(col2_label)}"], num_current_cols)
            write(file, join(sub_headers, " & ") * " \\\\\n")
            
            rules = []
            for i in 1:num_current_cols
                start_col = 3 + (i - 1) * 2
                end_col = start_col + 1
                push!(rules, "\\cmidrule(lr){$(start_col)-$(end_col)}")
            end
            write(file, join(rules, " ") * "\n")


            for m_idx in 1:n_rows
                current_m = m_start + m_idx - 1
                
                write(file, "$(row_labels[1]) & \\multirow{3}{*}{$(current_m)} & ")
                data_row_1 = []
                for k_idx in current_k_indices
                    stats = mat[m_idx, k_idx][1] # low_stats
                    push!(data_row_1, _format_latex_value(stats, column1, decimal_places))
                    push!(data_row_1, _format_latex_value(stats, column2, decimal_places))
                end
                write(file, join(data_row_1, " & ") * " \\\\\n")

                write(file, "$(row_labels[2]) & & ") # M-column is blank
                data_row_2 = []
                for k_idx in current_k_indices
                    stats = mat[m_idx, k_idx][2] # no_delay_stats
                    push!(data_row_2, _format_latex_value(stats, column1, decimal_places))
                    push!(data_row_2, _format_latex_value(stats, column2, decimal_places))
                end
                write(file, join(data_row_2, " & ") * " \\\\\n")

                write(file, "$(row_labels[3]) & & ") # M-column is blank
                data_row_3 = []
                for k_idx in current_k_indices
                    stats = mat[m_idx, k_idx][3] # high_stats
                    push!(data_row_3, _format_latex_value(stats, column1, decimal_places))
                    push!(data_row_3, _format_latex_value(stats, column2, decimal_places))
                end
                write(file, join(data_row_3, " & ") * " \\\\\n")

                if m_idx < n_rows
                    write(file, raw"\midrule" * "\n")
                end
            end

            write(file, raw"\bottomrule" * "\n")
            write(file, raw"\end{tabular}" * "\n")

            if block_idx < num_blocks
                write(file, raw"\vspace{1em}" * "\n\n") # 1em vertical space
            end
        end

        write(file, raw"\end{table}" * "\n")
    end # close file
    
    println("Successfully generated full dual-column LaTeX table at $(filename)")
end

function simple_9_node_network_stats()
    Random.seed!(SEED)
    simple = SDNUtils.SdnGraphUtils.simple_network(true)
    row_labels = ("BCC=270, BSC=350", "No Delay", "BCC=450, BSC=250")

    stats = mixed_strategies_stats(simple.g, 270.0, 450.0, 2000.0, 5, 4; k_low=2)

    stats_to_latex(stats, "results/simple_graph_payoff.tex", "Payoffs V* for the Simple 9-node network", :payoff; k_start=2, row_labels = row_labels)

    stats_to_latex_dual(stats, "results/simple_graph_placements.tex", "Number of non-zero probability controller placements and attack placements for the Simple 9-node network", :n_placements, "Placements", :n_attacks, "Attacks"; k_start=2, row_labels=row_labels)

    stats_to_latex_dual(stats, "results/simple_graph_exec_time.tex", "Execution time of the column-generation process for the 9-node network", :exec_time, "seconds", :iterations, "iterations"; k_start=2, row_labels=row_labels)

    stats_to_latex_dual(stats, "results/simple_graph_entropy.tex", "Entropy (bits) for the 9-node network", :placement_entropy, "Placements", :attack_entropy, "Attacks"; k_start=2, row_labels=row_labels)

    simple, stats
end

function cost266_network_stats()
    Random.seed!(SEED)
    cost266 = SDNUtils.read_sndgraph("graphs/cost266.sndlib"; weighted=true)
    row_labels = ("BCC=1500, BSC=1529.28", "No Delay", "BCC=2000, BCC=2000")

    stats = mixed_strategies_stats(cost266.g, 1500.0, 2000.0, 2000.0, 15, 6; k_low=2)

    stats_to_latex(stats, "results/cost266_graph_payoff.tex", "Payoffs V* for the cost266", :payoff; k_start=2, row_labels = row_labels)

    stats_to_latex_dual(stats, "results/cost266_graph_placements.tex", "Number of non-zero probability controller placements and attack placements for the cost266", :n_placements, "Placements", :n_attacks, "Attacks"; k_start=2, row_labels=row_labels)

    stats_to_latex_dual(stats, "results/cost266_graph_exec_time.tex", "Execution time of the column-generation process for the cost266", :exec_time, "seconds", :iterations, "iterations"; k_start=2, row_labels=row_labels)

    stats_to_latex_dual(stats, "results/cost266_graph_entropy.tex", "Entropy (bits) for the cost266", :placement_entropy, "Placements", :attack_entropy, "Attacks"; k_start=2, row_labels=row_labels)

    cost266, stats
end
