# vertices in the graph
set VERTICES;
# attacks. These are simply indices
set ATTACKS;
# Nodes affected by the attack
set VERTICES_A {ATTACKS} within VERTICES;
# vertices exceeding BCC {{v, w} in V^|2| : d(v, w) > BCC}
set U within {VERTICES, VERTICES};
# C(a) is the set of components resulting from an attack a
## Component Ids hold the indices of components resulting from attack a
set COMPONENT_IDS {ATTACKS};
set C {a in ATTACKS, COMPONENT_IDS[a]} within VERTICES;
# Sets related to delays
set W {VERTICES} within VERTICES;


# Which node is a primary controller
var y {VERTICES} binary;
# Which node is a backup controller
var x {VERTICES} binary;
# Component S survives attack a or not
var S {a in ATTACKS, COMPONENT_IDS[a]} binary;

var Y {ATTACKS} >= 0;

# Total number of primary controller nodes
param P_low >= 1 integer;
param P_high >= 1 integer;
# Total number of backup controller nodes
param B_low >= 0 integer;
param B_high >= 0 integer;

# Attack probability
param p_star {ATTACKS};

# (1a) Operator payoff
maximize OperatorPayoff:
    sum {a in ATTACKS} (p_star[a] * Y[a]);

# (1b)
s.t. NumberOfPrimaryControllers:
    P_low <= sum {v in VERTICES} y[v] <= P_high;

# (1c)
s.t. NumberOfBackupControllers:
    B_low <= sum {v in VERTICES} x[v] <= B_high;

# (1d) Respect Controller to controller delay for primary controllers only
s.t. SatisfyCCDelay {(v, w) in U}:
    y[v] + y[w] <= 1;

# (1e) Respect Switch to controller delay for primary controllers only
s.t. SatisfySCDelay {v in VERTICES}:
    sum {w in W[v]} y[w] >= 1;

# (1f) A node can only be one controller type
s.t. EnforceSingleControllerType {v in VERTICES}:
    y[v] + x[v] <= 1;

# (1g) Zero the components that  do not have a controller
s.t. ZeroComponentsWithoutControllers {a in ATTACKS, c_id in COMPONENT_IDS[a]}:
    S[a, c_id] <= sum {v in C[a, c_id]} (y[v] + x[v]);

# (1h) Count the number of surviving nodes
s.t. NumberOfSurvivingNodes {a in ATTACKS}:
    Y[a] = sum {c_id in COMPONENT_IDS[a]} card(C[a, c_id]) * S[a, c_id];
