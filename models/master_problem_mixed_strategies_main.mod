set A; # The set of attacks. Essentially, these are just indices.
set S; # The set of placements.

# V(s, a)
param V{S, A}; # How many nodes survive the attack

var y;
# q_s
var q{S} >= 0;

maximize Payoff:
    y;

s.t. EnsureProbabilityDistribution:
    sum {s in S} q[s] = 1.0;

s.t. BoundPayoff {a in A}:
    y <= sum {s in S} (V[s, a] * q[s]);
