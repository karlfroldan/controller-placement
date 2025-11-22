
# vertices in the graph
set VERTICES;
# attacks. These are simply indices
set ATTACKS;
# Nodes affected by the attack
set V_A {ATTACKS} within VERTICES;
# vertices exceeding BCC {{v, w} in V^|2| : d(v, w) > BCC}
set U within {VERTICES, VERTICES};
# C(a) is the set of components resulting from an attack a
## Component Ids hold the indices of components resulting from attack a
set COMPONENT_IDS {ATTACKS};
set C_A {a in ATTACKS, COMPONENT_IDS[a]} within VERTICES;

# Number of surviving nodes given attack
var Y {ATTACKS} integer >= 0;
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

s.t. ZeroComponentsWithoutControllers {a in ATTACKS, c_id in COMPONENT_IDS[a]}:
    sum {v in C_A[a, c_id]} y[v, a] <= card(C_A[a, c_id]) * sum {v in C_A[a, c_id]} s[v];

s.t. NumberOfSurvivingNodes {a in ATTACKS}:
    Y[a] = sum {v in VERTICES} y[v, a];
