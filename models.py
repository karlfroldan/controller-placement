import numpy as np

from itertools import product
from amplpy import AMPL

from utils import normalize_probabilities, set_range, values_to_list, to_set_ids

def component_ids_on_attack(network, attacks):
    # Contains the IDS of the components
    component_ids = {}
    # Components are indexed by the attack ID and the component ID after the attack.
    C_A = {}
    for i, a in attacks.items():
        cs = network.attack(a)
        component_ids[i] = set_range(cs)
        for j, c in enumerate(cs):
            C_A[(i, j)] = set(c)
    return component_ids, C_A

def pairs_exceeding_bcc_delay(network, bcc):
    return {
        (v, w) for v, w in product(network.nodes, network.nodes)
        if v < w and network.delays[v, w] > bcc
    }

def nodes_within_bsc_delay(network, bsc):
    W_V = {}
    for v in network.nodes:
        within_delay_list = {
            w for w in network.nodes if network.delays[v, w] <= bsc
        }
        W_V[v] = within_delay_list
    return W_V

class MathematicalModel:
    def __init__(self, model_file):

        assert(model_file)
        self.ampl = AMPL()
        # self.ampl.setOption('cplex_options', 'mipdisplay=0 lpmethod=0 barrier display=0 simplex display=0')
        # self.ampl.setOption('solver_msg', 0)

        # self.ampl.setOption('solver_msg', 0)
        # self.ampl.setOption('show_presolve_messages', 0)
        # self.ampl.setOption('cplex_options', 'mipdisplay=0 outlev=0 version=0')

        self._solved = False
        self._model_file = model_file
        self._model_loaded = False

        self._load_model()

    def _reset(self):
        self.ampl.eval('reset data;')

    def _load_model(self, model_file = None):
        if not self._model_loaded:
            f = model_file if model_file else self._model_file
            self.ampl.read(f)
            self._model_loaded = True

    def load_data(self):
        self._reset()
    
    def solve(self, solver='cplex', model_file = None):
        assert(self._model_loaded)

        self.ampl.option['solver'] = solver
        self.ampl.solve()
        self._solved = True

    def __str__(self):
        vars_map = self.ampl.get_variables()
        cons_map = self.ampl.get_constraints()

        if self._name:
            return f'Model {self._name}: {vars_map.size()} variables and {cons_map.size()} constraints'
        else:
            return 'Unnamed Model {vars_map.size()} variables and {cons_map.size()} constraints'

class MixedStrategyMasterProblem(MathematicalModel):
    def __init__(self, network):
        self._model_file = 'models/master_problem_mixed_strategies_main.mod';

        super(MixedStrategyMasterProblem, self).__init__(self._model_file)
        self._name = 'MixedStrategyMasterProblem'
        self.network = network

    def report(self):
        if not self._solved:
            self.solve()

        # This is the probability distribution of the attacker mixed strategy
        ps = self.ampl.get_constraint('BoundPayoff').get_values()
        ps = values_to_list(ps)

        qs = self.ampl.get_variables()['q']
        # Return the controller placement mixed-strategy.
        qs = values_to_list(qs)

        # payoff
        y = self.ampl.get_variable('y').value()

        x_constraint = self.ampl.get_constraint('EnsureProbabilityDistribution')
        x = x_constraint.dual()
        return {
            'y': y,
            'x': x,
            'q': normalize_probabilities(qs),
            'p': normalize_probabilities(ps),
        }

    def load_data(self, placements, attacks):
        super(MixedStrategyMasterProblem, self).load_data()
        # Convert them to IDs so we can pass them as sets.
        placement_mapping = to_set_ids(placements)
        attack_mapping = to_set_ids(attacks)

        self.ampl.set['S'] = sorted(placement_mapping.keys())
        self.ampl.set['A'] = sorted(attack_mapping.keys())

        v_data = {}

        for i, p in placement_mapping.items():
            for j, a in attack_mapping.items():
                surviving_nodes = len(self.network.surviving_nodes(p, a))
                assert(isinstance(surviving_nodes, int))

                v_data[(i, j)] = surviving_nodes

        self.ampl.param['V'] = v_data

class AttackGenerationPricingProblem(MathematicalModel):
    def __init__(self, network, K, eps = 1e-9):
        self._model_file = 'models/attack_generation_pricing_problem.mod'
        super(AttackGenerationPricingProblem, self).__init__(self._model_file)
        self.network = network

        self._name = 'AttackGenerationPricingProblem'

        self.eps = eps
        self.K = K

    def non_zero_placements(self, placements, q_star):
        new_placements = []
        new_q_star = []

        # Do not consider placements that have 0 probability.
        for p, q in zip(placements, q_star):
            if q >= self.eps:
                new_placements.append(p)
                new_q_star.append(q)
        return new_placements, new_q_star

    def report(self):
        if not self._solved:
            self.solve()

        var_a = self.ampl.get_variables()['a']
        a_star = values_to_list(var_a)
        return {
            'a*': a_star,
        }
        

    def load_data(self, placements, q_star):
        super(AttackGenerationPricingProblem, self).load_data()

        placements, q_star = self.non_zero_placements(placements, q_star)
        # Convert them to IDs
        s_prime = placements
        # s_prime = [i for i in range(len(self.placements))]
        placement_mapping = to_set_ids(s_prime)
        vertex_set = set(self.network.nodes)
        # vertex_set = list(self.network.nodes())
        edge_set = set(self.network.edges)

        v_s = placement_mapping

        self.ampl.set['V'] = vertex_set
        self.ampl.set['S_PRIME'] = set(placement_mapping.keys())
        self.ampl.set['E'] = edge_set
        self.ampl.set['V_S'] = v_s

        self.ampl.param['K'] = self.K
        self.ampl.param['q_star'] = q_star

class ControllerPlacementPricingProblem(MathematicalModel):
    def __init__(self, network, M, eps=1e-9):
        self._model_file = 'models/controller_placement_pricing_problem.mod'
        super(ControllerPlacementPricingProblem, self).__init__(self._model_file)

        self._name = 'ControllerPlacementPricingProblem'
        self.network = network
        self.M = M
        self.eps = eps

    def report(self):
        if not self._solved:
            self.solve()
        var_s = self.ampl.get_variables()['s']
        s_star = values_to_list(var_s)
        return {
            's*': s_star,
        }

    def load_data(self, attacks, p_star):
        super(ControllerPlacementPricingProblem, self).load_data()
        vertex_list = set(self.network.nodes)

        # Key => attack set pair
        attack_dict = to_set_ids(attacks)

        self.ampl.set['VERTICES'] = vertex_list
        self.ampl.set['ATTACKS'] = attack_dict.keys()
        
        V_A = {}
        for i, a in attack_dict.items():
            V_A[i] = a
        self.ampl.set['V_A'] = V_A

        C_IDS, C_A = component_ids_on_attack(self.network, attack_dict)

        self.ampl.set['COMPONENT_IDS'] = C_IDS
        self.ampl.set['C_A'] = C_A

        self.ampl.param['M'] = self.M
        self.ampl.param['p_star'] = p_star

class ControllerPlacementPricingProblemWithDelay(MathematicalModel):
    def __init__(self, network, M, bsc, bcc, eps=1e-9):
        self._model_file = 'models/controller_placement_pricing_problem_with_delay.mod'
        super(ControllerPlacementPricingProblemWithDelay, self).__init__(self._model_file)

        self._name = 'ControllerPlacementPricingProblem with Delay'
        self.network = network
        self.M = M
        self.eps = eps
        # Bound on the BSC delay
        self.bsc = bsc
        # Bound on the BCC delay
        self.bcc = bcc

    def report(self):
        if not self._solved:
            self.solve()
        var_s = self.ampl.get_variables()['s']
        s_star = values_to_list(var_s)
        # Y_star = values_to_list(self.ampl.get_variables()['Y'])

        V_star = self.ampl.get_objectives()['OperatorPayoff'].value()
        return {
            's*': s_star,
        }

    def load_data(self, attacks, p_star):
        super(ControllerPlacementPricingProblemWithDelay, self).load_data()
        vertex_list = set(self.network.nodes)

        # Key => attack set pair
        attack_dict = to_set_ids(attacks)

        self.ampl.set['VERTICES'] = vertex_list
        self.ampl.set['ATTACKS'] = attack_dict.keys()

        V_A = {}
        for i, a in attack_dict.items():
            # vs = self.network.remaining_nodes(a)
            V_A[i] = a
        self.ampl.set['V_A'] = V_A

        self.ampl.set['U'] = pairs_exceeding_bcc_delay(self.network, self.bcc)
        self.ampl.set['W_V'] = nodes_within_bsc_delay(self.network, self.bsc)
        
        # Components after the attack. c in C(a) contains the list of vertices in each component
        C_IDS, C_A = component_ids_on_attack(self.network, attack_dict)

        self.ampl.set['COMPONENT_IDS'] = C_IDS
        self.ampl.set['C_A'] = C_A

        self.ampl.param['M'] = self.M
        self.ampl.param['p_star'] = p_star

# PURE IMPLEMENTATIONS
class CPOP(MathematicalModel):
    def __init__(self, network):
        self._model_file = 'models/cpop.mod'
        super(CPOP, self).__init__(self._model_file)
        self._name = 'CPOP'

        self.network = network

    def report(self):
        if not self._solved:
            self.solve()

        var_s = self.ampl.get_variables()['s']
        objective_Y = self.ampl.get_objective('SurvivingNodes')
        s_star = values_to_list(var_s)
        return {
            's': s_star,
            'Y': objective_Y.value(),
        }

    def load_data(self, attacks, M):
        super(CPOP, self).load_data()
        vertex_list = set(self.network.nodes)

        # Key => attack set pair
        attack_dict = to_set_ids(attacks)

        V_A = {}
        for i, a in attack_dict.items():
            # vs = self.network.remaining_nodes(a)
            V_A[i] = a

        C_IDS, C_A = component_ids_on_attack(self.network, attack_dict)

        self.ampl.set['VERTICES'] = vertex_list
        self.ampl.set['ATTACKS'] = attack_dict.keys()
        self.ampl.set['COMPONENT_IDS'] = C_IDS
        self.ampl.set['C_A'] = C_A
        self.ampl.set['V_A'] = V_A
        
        self.ampl.param['M'] = M

class NAOP(MathematicalModel):
    def __init__(self, network):
        self.network = network
        self._model_file = 'models/naop.mod'
        super(NAOP, self).__init__(self._model_file)
        self.name = 'NAOP'

    def load_data(self, placements, K):
        super(NAOP, self).load_data()

        placement_ids = set_range(placements)

        vertex_list = set(self.network.nodes)

        edge_ids = set_range(self.network.edges)

        edgeid2edge = {i: e for i, e in enumerate(self.network.edges)}
        edge2edgeid = {}

        for i, e in enumerate(self.network.edges):
            edge2edgeid[e] = i
            edge2edgeid[(e[1], e[0])] = i

        # Set of links incident with node v
        delta = {}
        for v in vertex_list:
            incidents = self.network.edges(v)
            incident_ids = [edge2edgeid[e] for e in incidents]
            delta[v] = set(incident_ids)

        V_S = {i: set(s) for i, s in enumerate(placements)}

        alpha = {i: e[0] for i, e in enumerate(edgeid2edge)}
        beta = {i: e[1] for i, e in enumerate(edgeid2edge)}

        # SETS
        self.ampl.set['VERTICES'] = set(vertex_list)
        self.ampl.set['EDGES'] = set(edge_ids)
        self.ampl.set['DELTA'] = delta
        self.ampl.set['PLACEMENTS'] = set(placement_ids)
        self.ampl.set['V_S'] = V_S

        # PARAMS
        self.ampl.param['K'] = K
        self.ampl.param['alpha'] = alpha
        self.ampl.param['beta'] = beta

    def report(self):
        if not self._solved:
            self.solve()

        var_a = self.ampl.get_variables()['a']
        objective_Z = self.ampl.get_objective('SurvivingNodes')
        a_star = values_to_list(var_a)

        return {
            'a': a_star,
            'Z': objective_Z.value()
        }
