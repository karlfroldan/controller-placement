### Create the most dangerous k-node topological attack

# Set of all vertices in the graph
set VERTICES;

set ATTACK_IDX;
set ATTACKS {ATTACK_IDX} within VERTICES;

# Directed superset of all vertices
set H within {VERTICES, VERTICES};

# The set of adjacent nodes for some node v
set ADJ {VERTICES} within VERTICES;

# The degree of a node in G
param deg {VERTICES} integer;

# \mathcal V(v, w) for (v, w) in H
set V_PAIR {(v, w) in H} =
    if deg[v] <= deg[w] then
        ADJ[v]
    else
        ADJ[w];

# The set of directed node-pairs representing the set of undirected link
# of graph G
set H_EDGES within H;
# Set complement E - H_EDGES
set H_PRIME_EDGES = H diff H_EDGES;

# How many nodes to attack
param K;

# The nodes to attack
var a{VERTICES} binary;
var u{H} binary;

minimize Payoff:
    sum {(v, w) in H} u[v, w];

s.t. KNodeAttackConstraint:
    sum {v in VERTICES} a[v] = K;

s.t. ConstraintA5C {(v, w) in H_EDGES}:
    u[v, w] + a[v] + a[w] >= 1;

# u_vt if v < t and u_tv if v > t - min(v, t), max(v, t)
s.t. ConstraintA5D {(v, w) in H_PRIME_EDGES, t in V_PAIR[v, w]}:
    u[v, w] >= u[min(v, t), max(v, t)] + u[min(t, w), max(t, w)] + a[t] - 1;

s.t. AttackUniquenessConstraint {aidx in ATTACK_IDX}:
    sum {v in ATTACKS[aidx]} a[v] <= K - 1;
