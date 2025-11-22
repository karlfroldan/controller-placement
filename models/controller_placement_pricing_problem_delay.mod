
# vertices in the graph
set VERTICES;
# attacks. These are simply indices
set ATTACKS;
# Nodes affected by the attack
set V_A {ATTACKS} within VERTICES;
# vertices exceeding BCC {{v, w} in V^|2| : d(v, w) > BCC}
set U within VERTICES cross VERTICES;
# V_C is the set of components resulting from an attack a
set V_C {ATTACKS} within VERTICES;
# The nodes within the SC-delay constraints of node v
set W_V {VERTICES} within VERTICES;

# Number of surviving nodes given attack
var Y {ATTACKS};
# Whether node v is a controller
var s {VERTICES} binary;
# The number of surviving nodes given attack a
var y {VERTICES, ATTACKS} binary;

# Total number of controller nodes
param M integer;
# Attack probability
param p_star {ATTACKS};

maximize OperatorPayoff:
    sum {a in ATTACKS} (p_star[a] * Y[a]);

s.t. NumberOfControllers:
    sum {v in VERTICES} s[v] = M;

s.t. ZeroAttackedNodes {a in ATTACKS, v in V_A[a]}:
    y[v, a] = 0;
