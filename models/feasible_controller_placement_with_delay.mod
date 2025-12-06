# Find a suitable controller placement with delays.
# This need not be optimal
set VERTICES;
set L;
# The controllers in placement l
set T {L} within VERTICES;

set W {VERTICES} within VERTICES;
set U within VERTICES cross VERTICES;

# How many controllers should there be
param big_M;
# Newly generated placement should differ from the already-generated one.
param small_m;

var Z >= 0;
var s {VERTICES} binary;

minimize Payoff:
    0;

s.t. NumberOfControllers:
    sum {v in VERTICES} s[v] = big_M;

s.t. SatisfyCCDelay {(v, w) in U}:
    s[v] + s[w] <= 1;

s.t. SatisfySCDelay {v in VERTICES}:
    sum {w in W[v]} s[w] >= 1;

s.t. NewPlacementMustDiffer {l in L}:
    sum {v in T[l]} s[v] <= big_M - small_m;
