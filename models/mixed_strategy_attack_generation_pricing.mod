# The number of nodes
set VERTICES;
# Mon-zero probability placements. Just a bunch of indices, actually.
set PLACEMENTS;
# The number of edges
set EDGES within VERTICES cross VERTICES;
# The nodes that are included in placement s
set V_S {PLACEMENTS} within VERTICES;

var F{PLACEMENTS} integer >= 0;
var a{VERTICES} binary;
var z{VERTICES, PLACEMENTS} binary;

# Payoff of the network operator
param q_star{PLACEMENTS};
# Total required attacked nodes.
param K_low;
param K_high;

param cost {VERTICES};
param budget;

minimize ControllerPayoff:
    sum {s in PLACEMENTS} (q_star[s] * F[s]);

s.t. ConstraintAttackSize:
    K_low <= sum {v in VERTICES} a[v] <= K_high;

s.t. ConstraintNodeSurvival {s in PLACEMENTS, v in V_S[s]}:
    z[v, s] = 1 - a[v];

s.t. ConstraintNeighborhoodSurvival1 {(alpha, beta) in EDGES, s in PLACEMENTS}:
    z[beta, s] >= z[alpha, s] - a[beta];

s.t. ConstraintNeighborhoodSurvival2 {(alpha, beta) in EDGES, s in PLACEMENTS}:
    z[alpha, s] >= z[beta, s] - a[alpha];

s.t. ConstraintNumberOfSurvivingNodes {s in PLACEMENTS}:
    F[s] = sum {v in VERTICES} z[v, s];

s.t. WithinCost:
    sum {v in VERTICES} (a[v] * cost[v]) <= budget;
