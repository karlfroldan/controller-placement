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
    
    m, D
end

function cpop(
	g :: SimpleGraph, 
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
	# First, let's get all the components resaulting from a given attack a
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
	g :: SimpleGraph,
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
    g :: SimpleGraph, M :: Int, K :: Int;
    TOL :: Float64 = 1e-6, silent = false
)

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

        time_begin = time()
		optimize!(naop_model)
        time_end = time()

		assert_is_solved_and_feasible(naop_model)

        push!(naop_times, time_end - time_begin)
		
		a_star = SDNUtils.to_placement(value.(a_star))
		Z_star = objective_value(naop_model)

		if Z_star >= Y_star - TOL
			break
		end

        if !silent
            @printf("    [%d] POST-NAOP: Z* = %.2f, time = %.5f\n", n_iterations, Z_star,
                    time_end - time_begin)
        end

		# STEP 2: A = A union {a*}. Then solve P[M, A] to get placement s*. Repeat
		if !(a_star ∈ attack_set)
			push!(attack_set, a_star)
		else
            error("This should not happen");
		end

		cpop_model, s_star = cpop(g, M, attack_set)
        time_begin = time()
		optimize!(cpop_model)
        time_end = time()

		assert_is_solved_and_feasible(cpop_model)

		s_star = SDNUtils.to_placement(value.(s_star))
		Y_star = objective_value(cpop_model)
		n_iterations += 1

        if !silent
            @printf("    [%d] POST-CPOP: time = %.5f\n", n_iterations, time_end - time_begin)
            @printf("    [%d] a_star = %s\n", n_iterations, a_star)
            @printf("    [%d] s_star = %s\n", n_iterations, s_star)
            @printf("    [%d] Controller Placement Optimization: Y* = %.2lf, Z* = %.2lf\n", n_iterations, Y_star, Z_star)
        end

        push!(cpop_times, time_end - time_begin)
	end

    if !silent
        @printf("CONTROLLER PURE [%d] # of generated attacks: %d\n", n_iterations, length(attack_set))
        @printf("CONTROLLER PURE [%d] cpop time total, H-mean: (%.5f, %.5f), naop time total, H-mean: (%.5f, %.5f)\n",
                n_iterations, sum(cpop_times), harmmean(cpop_times),
                sum(naop_times), harmmean(naop_times))
    end

	s_star, Y_star
end

function attack_optimization(
    g :: SimpleGraph, M :: Int, K :: Int;
    TOL :: Float64 = 1e-6, silent = false
)

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
		cpop_model, s_star = cpop(g, M, attack_set)

		begin_time = time()
		optimize!(cpop_model)
		elapsed_time = time() - begin_time
		push!(cpop_times, elapsed_time)

		is_solved_and_feasible(cpop_model)

		s_star = SDNUtils.to_placement(value.(s_star))
		Y_star = objective_value(cpop_model)

		# println("[$n_iterations] Z_star = $Z_star, Y_star = $Y_star")
		if Y_star <= Z_star + TOL
			break
		end

		# STEP 2: S = S union {s*}. Then solve A[K, S] to get attack a*. Repeat
		if !(s_star ∈ placement_set)
			push!(placement_set, s_star)
		else
            error("This should not happen")
		end

		naop_model, a_star = naop(g, K, placement_set)
		begin_time = time()
		optimize!(naop_model)
		elapsed_time = time() - begin_time
		push!(naop_times, elapsed_time)

		is_solved_and_feasible(naop_model)

		a_star = SDNUtils.to_placement(value.(a_star))
		Z_star = objective_value(naop_model)

        if !silent
            @printf("    [%d] a_star = %s\n", n_iterations, a_star)
            @printf("    [%d] s_star = %s\n", n_iterations, s_star)
            @printf("    [%d] Attack Optimization: Y* = %.2lf, Z* = %.2lf\n", n_iterations, Y_star, Z_star)
        end

		n_iterations += 1
	end
    
    if !silent
        @printf("ATTACKER PURE [%d] cpop time total, H-mean: (%.5f, %.5f), naop time total, H-mean: (%.5f, %.5f)\n",
                n_iterations, sum(cpop_times), harmmean(cpop_times),
                sum(naop_times), harmmean(naop_times))
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

function mixed_strategy_algorithm(
	g :: AbstractGraph, 
	M :: Int, K :: Int;
    silent = false,
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
	s = SDNUtils.gen_random_vector(nv(g), M)
	a = SDNUtils.gen_random_vector(nv(g), K)

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

	# @show x_star, y_star

    placement_times = []
    attack_times = []
    n_iterations = 0

    has_updated = true

    if !silent
        @printf("RUNNING MIXED-STRATEGY FOR M=%d, K=%d\n", M, K)
    end

	while has_updated
        if !silent
            @printf("    ---- Iteration %d ----\n", n_iterations)
        end
        has_updated = false
		
		# Solve the placement generation problem to get a new placement s'
        if has_delays
            cp, s = pricing_placement(g, M, attacks, p, delay_table, bcc, bsc)
        else
		    cp, s = pricing_placement(g, M, attacks, p)
        end

        time_begin = time()
		optimize!(cp)
        time_end = time()

        push!(placement_times, time_end - time_begin)
        if !silent
            @printf("    [%d] pricing placement time: %.5f\n", n_iterations, time_end - time_begin)
        end
        
		assert_is_solved_and_feasible(cp)
		s = SDNUtils.to_placement(value.(s))
		# Generated a new placement s
		outcomes = [SDNUtils.game_outcome(g, s, a′) for a′ ∈ attacks]
		expected = outcomes ⋅ p

        if !silent
            @printf("    [%d] Expected Value = %.5f, x_star = %.5f\n", n_iterations, expected, x_star)
        end
		if expected > x_star - TOL
			if s ∉ placements
				has_updated = true
				push!(placements, s)
			end
		end

		# Then solve the master problem again.
		master, q, shadow_prices = master_placement(g, placements, attacks)
        time_begin = time()
		optimize!(master)
        time_end = time()

        if !silent
            @printf("    [%d] master time: %.5f\n", n_iterations, time_end - time_begin)
        end
        
		assert_is_solved_and_feasible(master)

		q = value.(q)
		p = SDNUtils.project_simplex(-dual.(shadow_prices))

		x_star = dual_objective_value(master)
		y_star = objective_value(master)

		# Then generate a new attack
		na, a = pricing_attack(g, K, placements, q)
        time_begin = time()
		optimize!(na)
        time_end = time()
        push!(attack_times, time_end - time_begin)

        if !silent
            @printf("    [%d] pricing attack time: %.5f\n", n_iterations, time_end - time_begin)
        end

		assert_is_solved_and_feasible(na)

		a = SDNUtils.to_placement(value.(a))
		outcomes = [SDNUtils.game_outcome(g, s′, a) for s′ ∈ placements]
		expected = outcomes ⋅ q

        if !silent
            @printf("    [%d] Expected Value = %.5f, y_star = %.5f\n", n_iterations, expected, y_star)
        end
		
		if y_star > expected - TOL
			if a ∉ attacks
				has_updated = true
				push!(attacks, a)
			end
		end

		master, q, shadow_prices = master_placement(g, placements, attacks)
        time_begin = time()
		optimize!(master)
        time_end = time()
		assert_is_solved_and_feasible(master)

        if !silent
            @printf("    [%d] master time: %.5f\n", n_iterations, time_end - time_begin)
        end

		q = value.(q)
		p = SDNUtils.project_simplex(-dual.(shadow_prices))
		x_star = dual_objective_value(master)
		y_star = objective_value(master)

        n_iterations += 1
	end

    if !silent
        @printf("MIXED-STRATEGIES [%d] placement time total, H-mean: (%.5f, %.5f), attack time total, H-mean: (%.5f, %.5f)\n",
                n_iterations, sum(placement_times), harmmean(placement_times),
                sum(attack_times), harmmean(attack_times))
    end

    (
        p_star = p,
        q_star = q,
        V_star = x_star,
        attacks = attacks,
        placements = placements
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

function outcome_metrics_pure(g, M::Int, K::Int)
    V = nv(g)
    @assert 0 ≤ M ≤ V && 0 ≤ K ≤ V
	_, v_maxmin = controller_placement_optimization(g, M, K)
	_, v_minmax = attack_optimization(g, M, K)

    return v_minmax, v_maxmin
end

function outcome_table(g, M_max :: Int, K_max :: Int; K_min :: Int = 2, full = false)
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
