set VERTICES;
# Set indices of edges
set EDGES;
# Set indices of placements
set PLACEMENTS;
set V_S {PLACEMENTS} within VERTICES;

# Incident vertices
set DELTA {VERTICES} within EDGES;

param K integer >= 0; # K-node attacks
param alpha {EDGES} in VERTICES;
param beta {EDGES} in VERTICES;

var Z >= 0;
# Whether a node v is attacked
var a {VERTICES} binary;
# t_e is not available after the attack
var t {EDGES} binary;
# Node v survives the attack on placement s
var z {VERTICES, PLACEMENTS} binary;

minimize SurvivingNodes:
    Z;

s.t. EnsureKNodeAttack:
    sum {v in VERTICES} a[v] = K;

s.t. AttackIncidentLinks {v in VERTICES, e in DELTA[v]}:
    t[e] >= a[v];

s.t. AttackIncidentNodes {e in EDGES}:
    t[e] <= a[alpha[e]] + a[beta[e]];

s.t. NodeSurvivalUpper {s in PLACEMENTS, v in VERTICES}:
    z[v, s] <= 1 - a[v];

s.t. NodeSurvivalLower {s in PLACEMENTS, v in V_S[s]}:
    z[v, s] >= 1 - a[v];

s.t. IncidentSurvival1 {s in PLACEMENTS, e in EDGES}:
    z[alpha[e], s] >= z[beta[e], s] - t[e];

s.t. IncidentSurvival2 {s in PLACEMENTS, e in EDGES}:
    z[beta[e], s] >= z[alpha[e], s] - t[e];

s.t. PayoffConstraint {s in PLACEMENTS}:
    Z >= sum {v in VERTICES} z[v, s];
