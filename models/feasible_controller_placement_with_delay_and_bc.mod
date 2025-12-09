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
param num_controllers >= 1 integer;
param P_low >= 1 integer;
param P_high >= 1 integer;
param B_low >= 0 integer;
param B_high >= 0 integer;
# Newly generated placement should differ from the already-generated one.
param m integer;

# Big-M
param M_primary := m;
param M_backup := m;

var Z >= 0;

var y {VERTICES} binary;
var x {VERTICES} binary;

# Conditional variables
var y_cond_primary binary;
var y_cond_backup binary;

minimize Payoff:
    0;

s.t. NumberOfTotalControllers:
    sum {v in VERTICES} (y[v] + x[v]) <= num_controllers;

s.t. NumberOfControllersPrimary:
    P_low <= sum {v in VERTICES} y[v] <= P_high;

s.t. NumberOfControllersBackup:
    B_low <= sum {v in VERTICES} x[v] <= B_high;

s.t. SatisfyCCDelay {(v, w) in U}:
    y[v] + y[w] <= 1;

s.t. SatisfySCDelay {v in VERTICES}:
    sum {w in W[v]} y[w] >= 1;

s.t. EnforceSingleNodeType {v in VERTICES}:
    y[v] + x[v] <= 1;

# We want to enforce at least one of these to be true.
s.t. NewPlacementMustDifferPrimary {l in L}:
    sum {v in T_primary[l]} y[v] <= (sum {v in VERTICES} y[v]) - m + (M_primary * y_cond_primary);

s.t. NewPlacementMustDifferBackup {l in L}:
    sum {v in T_backup[l]} x[v] <= (sum {v in VERTICES} x[v]) - m + (M_backup * y_cond_backup);

s.t. EnforceAtLeastOneSatisfies:
    y_cond_primary + y_cond_backup <= 1;
