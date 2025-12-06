# Find a suitable controller placement with delays.
# This need not be optimal
set VERTICES;
set L;
# The controllers in placement l
set T_primary {L} within VERTICES;
set T_backup {L} within VERTICES;

set W {VERTICES} within VERTICES;
set U within VERTICES cross VERTICES;

# How many controllers should there be
param P;
param B;
# Newly generated placement should differ from the already-generated one.
param m;

# Big-M
param M_primary := m;
param M_backup := m;

var Z >= 0;


var y {VERTICES} binary;
var x {VERTICES} binary;

# Conditional/Relaxation variables
var y_cond_primary binary;
var y_cond_backup binary;

minimize Payoff:
    0;

s.t. NumberOfControllersPrimary:
    sum {v in VERTICES} y[v] = P;

s.t. NumberOfControllersBackup:
    sum {v in VERTICES} x[v] = B;

s.t. SatisfyCCDelay {(v, w) in U}:
    y[v] + y[w] <= 1;

s.t. SatisfySCDelay {v in VERTICES}:
    sum {w in W[v]} y[w] >= 1;

s.t. EnforceSingleNodeType {v in VERTICES}:
    y[v] + x[v] <= 1;

# We want to enforce at least one of these to be true.
s.t. NewPlacementMustDifferPrimary {l in L}:
    sum {v in T_primary[l]} y[v] <= (P - m) + (M_primary * y_cond_primary);

s.t. NewPlacementMustDifferBackup {l in L}:
    sum {v in T_backup[l]} x[v] <= (B - m) + (M_backup + y_cond_backup);

s.t. EnforceAtLeastOneSatisfies:
    y_cond_primary + y_cond_backup <= 1;
