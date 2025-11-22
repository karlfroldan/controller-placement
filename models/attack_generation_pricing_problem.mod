set V; # The number of nodes
set S_PRIME; # Mon-zero probability placements. Just a bunch of indices, actually.
set E within {V, V}; # The number of edges
set V_S {S_PRIME} within V; # Controller nodes for s in S_PRIME

var F{S_PRIME} integer >= 0;
var a{V} binary;
var z{S_PRIME, V} binary;

# Payoff of the network operator
param q_star{S_PRIME};
# Total required attacked nodes.
param K;

minimize ControllerPayoff:
    sum {s in S_PRIME} (q_star[s] * F[s]);

s.t. ConstraintAttackSize:
    sum {v in V} a[v] = K;

s.t. ConstraintNodeSurvival {s in S_PRIME, v in V_S[s]}:
    z[s, v] = 1 - a[v];

s.t. ConstraintNeighborhoodSurvival1 {(alpha, beta) in E, s in S_PRIME}:
    z[s, beta] >= z[s, alpha] - a[beta];

s.t. ConstraintNeighborhoodSurvival2 {(alpha, beta) in E, s in S_PRIME}:
    z[s, alpha] >= z[s, beta] - a[alpha];

s.t. ConstraintNumberOfSurvivingNodes {s in S_PRIME}:
    F[s] = sum {v in V} z[s, v];
