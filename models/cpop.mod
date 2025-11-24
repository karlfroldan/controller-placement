
# Set of attacks
set ATTACKS;
set VERTICES;
# Set of nodes to attack given attack a
set V_A {ATTACKS} within VERTICES;
# C(a) is the set of components resulting from an attack a
## Component Ids hold the indices of components resulting from attack a
set COMPONENT_IDS {ATTACKS};
set C_A {a in ATTACKS, COMPONENT_IDS[a]} within VERTICES;

# How many controllers to use
param M;

var Y >= 0;
var s {VERTICES} binary;
var y {VERTICES, ATTACKS} binary;

maximize SurvivingNodes:
    Y;

s.t. NumberOfControllers:
    sum {v in VERTICES} s[v] = M;

s.t. ZeroAttackedNodes {a in ATTACKS, v in V_A[a]}:
    y[v, a] = 0;

s.t. ZeroComponentsWithoutControllers {a in ATTACKS, c_id in COMPONENT_IDS[a]}:
    sum {v in C_A[a, c_id]} y[v, a] <= card(C_A[a, c_id]) * sum {v in C_A[a, c_id]} s[v];

s.t. NumberOfSurvivingNodes {a in ATTACKS}:
    Y <= sum {v in VERTICES} y[v, a];
