
set V; # The number of nodes
set S_PRIME; # Mon-zero probability placements
set E; # The number of edges

var F{S_PRIME} integer >= 1;
var a{V} binary;
var z{S_PRIME, V} binary;

# Payoff of the network operator
param q_star{S_PRIME};
# Total required attacked nodes.
param K;

min ControllerPayoff:
    sum {s in S_PRIME} (q_star[s] * F[s]);

s.t. ConstraintAttackSize:
    sum {v in V} a[v] = K;

s.t. ConstraintNodeSurvival {v in V, s in S_PRIME}:
    z[s, v] = 1 - a[v];

# s.t. ConstraintNeighborhoodSurvival {e in EDGES
