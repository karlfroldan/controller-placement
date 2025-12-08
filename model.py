import numpy as np

from itertools import product
from amplpy import AMPL

import networkx as nx

from utils import normalize_probabilities, set_range, values_to_list, to_set_ids, one_indices

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

        self._solved = False
        self._model_file = model_file
        self._model_loaded = False

        self._load_model()

    def _reset(self):
        self._solved = False
        self.ampl.eval('reset data;')

    def _load_model(self, model_file = None):
        if not self._model_loaded:
            f = model_file if model_file else self._model_file
            self.ampl.read(f)
            self._model_loaded = True

    def get_result_array(self, v):
        vs = self.ampl.get_variables()[v]
        vs = values_to_list(vs)
        return vs

    def load_data(self):
        self._reset()
    
    def solve(self, solver='cplex', model_file = None):
        assert(self._model_loaded)

        self.ampl.option['solver'] = solver
        self.ampl.solve()
        self._solved = True

        # print(f'Result = {self.ampl.solve_result}')
        assert self.ampl.solve_result == "solved", self.ampl.solve_result


    def __str__(self):
        vars_map = self.ampl.get_variables()
        cons_map = self.ampl.get_constraints()

        if self._name:
            return f'Model {self._name}: {vars_map.size()} variables and {cons_map.size()} constraints'
        else:
            return 'Unnamed Model {vars_map.size()} variables and {cons_map.size()} constraints'

class MixedStrategyMasterProblem(MathematicalModel):
    def __init__(self, network):
        self._model_file = 'models/mixed_strategy_master_problem.mod';

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
        # We have backup controllers
        if isinstance(placements, tuple):
            primary_controllers = placements[0]
            backup_controllers = placements[1]

            primary_placement_mapping = to_set_ids(primary_controllers)
            backup_placement_mapping = to_set_ids(backup_controllers)

            placement_keys = primary_placement_mapping.keys()
        else:
            placement_mapping = to_set_ids(placements)
            placement_keys = placement_mapping.keys()

        # Convert them to IDs so we can pass them as sets.

        attack_mapping = to_set_ids(attacks)

        self.ampl.set['S'] = sorted(placement_keys)
        self.ampl.set['A'] = sorted(attack_mapping.keys())

        v_data = {}

        for i in placement_keys:
            for j, a in attack_mapping.items():
                if isinstance(placements, tuple):
                    p = primary_placement_mapping[i]
                    b = backup_placement_mapping[i]
                    surviving_nodes = self.network.surviving_nodes(p, a, backup_controllers = b)
                else:
                    surviving_nodes = self.network.surviving_nodes(p, a)
                n_survivors = len(surviving_nodes)
                
                assert(isinstance(n_survivors, int))

                v_data[(i, j)] = n_survivors

        self.ampl.param['V'] = v_data

class AttackGenerationPricingProblem(MathematicalModel):
    def __init__(self, network, K, eps = 1e-9, budget=0, costs=None):
        self._model_file = 'models/mixed_strategy_attack_generation_pricing.mod'
        super(AttackGenerationPricingProblem, self).__init__(self._model_file)
        self.network = network

        self._name = 'AttackGenerationPricingProblem'

        self.eps = eps
        self.K = K
        self.budget = budget
        self.costs = costs

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
        """
        Possible values for `costs` can be `degree`, `centrality`, or a table of costs for each node.
        If set as None, then all costs will be 0.
        """
        super(AttackGenerationPricingProblem, self).load_data()

        # Has backup controllers
        if isinstance(placements, tuple):
            # Combine the primary and backup controllers
            
            placements = [sorted(pc + bc) for pc, bc in zip(placements[0], placements[1])]

        placements, q_star = self.non_zero_placements(placements, q_star)

        # Convert them to IDs
        s_prime = placements
        placement_mapping = to_set_ids(s_prime)
        # print(placement_mapping)
        vertex_set = set(self.network.nodes)
        edge_set = set(self.network.edges)

        self.ampl.set['VERTICES'] = vertex_set
        self.ampl.set['PLACEMENTS'] = set(placement_mapping.keys())
        self.ampl.set['EDGES'] = edge_set
        # How placement IDs related to the actual placements.
        self.ampl.set['V_S'] = placement_mapping

        if isinstance(self.K, int):
            self.ampl.param['K_low'] = self.K
            self.ampl.param['K_high'] = self.K
        else:
            self.ampl.param['K_low'] = self.K[0]
            self.ampl.param['K_high'] = self.K[1]
        self.ampl.param['q_star'] = q_star

        self.ampl.param['budget'] = self.budget
        
        if self.costs is None:
            self.ampl.param['cost'] = to_set_ids([0 for _ in vertex_set])
        elif self.costs == 'degree':
            assert(self.budget > 0)
            cost_table = { v: self.network.degree(v) for v in vertex_set }
            self.ampl.param['cost'] = cost_table
            # print(f'cost_table: {cost_table}')
        elif self.costs == 'load_centrality':
            assert(self.budget > 0)
            cost_table = nx.load_centrality(self.network.g)
            # cost_table = { k: v for k, v in nx.load_centrality(self.network.g).items() }
            self.ampl.param['cost'] = cost_table
            # print(f'cost_table: {cost_table}')
        else:
            assert(self.budget > 0)
            self.ampl.param['cost'] = cost_table

class ControllerPlacementPricingProblem(MathematicalModel):
    def __init__(self, network, M, eps=1e-9):
        self._model_file = 'models/mixed_strategy_controller_placement_pricing.mod'
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


        
        V_A = {}
        for i, a in attack_dict.items():
            V_A[i] = a

        C_IDS, C_A = component_ids_on_attack(self.network, attack_dict)

        self.ampl.set['VERTICES'] = vertex_list
        self.ampl.set['ATTACKS'] = attack_dict.keys()
        self.ampl.set['COMPONENT_IDS'] = C_IDS
        self.ampl.set['C_A'] = C_A
        self.ampl.set['V_A'] = V_A

        self.ampl.param['M'] = self.M
        self.ampl.param['p_star'] = p_star

class ControllerPlacementPricingProblemWithDelay(MathematicalModel):
    def __init__(self, network, M, bsc, bcc, eps=1e-9):
        self._model_file = 'models/mixed_strategy_controller_placement_pricing_with_delay.mod'
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
        V_A = {}
        for i, a in attack_dict.items():
            # vs = self.network.remaining_nodes(a)
            V_A[i] = a
        
        # Components after the attack. c in C(a) contains the list of vertices in each component
        C_IDS, C_A = component_ids_on_attack(self.network, attack_dict)

        self.ampl.set['VERTICES'] = vertex_list
        self.ampl.set['ATTACKS'] = attack_dict.keys()
        self.ampl.set['U'] = pairs_exceeding_bcc_delay(self.network, self.bcc)
        self.ampl.set['VERTICES_A'] = V_A
        self.ampl.set['W'] = nodes_within_bsc_delay(self.network, self.bsc)
        self.ampl.set['COMPONENT_IDS'] = C_IDS
        self.ampl.set['C'] = C_A

        self.ampl.param['M'] = self.M
        self.ampl.param['p_star'] = p_star

class ControllerPlacementPricingProblemWithDelayAndBC(MathematicalModel):
    def __init__(self, network, P, B, bsc, bcc, eps=1e-9):
        self._model_file = 'models/mixed_strategy_controller_placement_pricing_with_delay_bc.mod'
        super().__init__(self._model_file)

        self._name = 'ControllerPlacementPricingProblem with Delay and Backup Controllers'
        self.network = network
        self.eps = eps
        # Bound on the BSC delay
        self.bsc = bsc
        # Bound on the BCC delay
        self.bcc = bcc

        self.P = P
        self.B = B

    def report(self):
        if not self._solved:
            self.solve()
        var_primary = self.ampl.get_variables()['y']
        y_star = values_to_list(var_primary)
        # Y_star = values_to_list(self.ampl.get_variables()['Y'])

        # V_star = self.ampl.get_objectives()['OperatorPayoff'].value()

        var_backup = self.ampl.get_variables()['x']
        x_star = values_to_list(var_backup)
        
        return {
            'y*': y_star,
            'x*': x_star,
        }

    def load_data(self, attacks, p_star):
        super().load_data()
        vertex_list = set(self.network.nodes)

        # Key => attack set pair
        attack_dict = to_set_ids(attacks)
        V_A = {}
        for i, a in attack_dict.items():
            # vs = self.network.remaining_nodes(a)
            V_A[i] = a
        
        # Components after the attack. c in C(a) contains the list of vertices in each component
        C_IDS, C_A = component_ids_on_attack(self.network, attack_dict)

        self.ampl.set['VERTICES'] = vertex_list
        self.ampl.set['ATTACKS'] = attack_dict.keys()
        self.ampl.set['U'] = pairs_exceeding_bcc_delay(self.network, self.bcc)
        self.ampl.set['VERTICES_A'] = V_A
        self.ampl.set['W'] = nodes_within_bsc_delay(self.network, self.bsc)
        self.ampl.set['COMPONENT_IDS'] = C_IDS
        self.ampl.set['C'] = C_A

        if isinstance(self.P, int):
            self.ampl.param['P_low'] = self.P
            self.ampl.param['P_high'] = self.P
        else:
            self.ampl.param['P_low'] = self.P[0]
            self.ampl.param['P_high'] = self.P[1]

        if isinstance(self.B, int):
            self.ampl.param['B_low'] = self.B
            self.ampl.param['B_high'] = self.B
        else:
            self.ampl.param['B_low'] = self.B[0]
            self.ampl.param['B_high'] = self.B[1]            
        self.ampl.param['p_star'] = p_star

# PURE IMPLEMENTATIONS
class CPOP(MathematicalModel):
    def __init__(self, network):
        self._model_file = 'models/pure_strategy_cpop.mod'
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
        self._model_file = 'models/pure_strategy_naop.mod'
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

        # print(edgeid2edge)
        alpha = {i: e[0] for i, e in edgeid2edge.items()}
        beta = {i: e[1] for i, e in edgeid2edge.items()}

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

class FeasibleControllerPlacementWithDelay(MathematicalModel):
    def __init__(self, network, M, bcc, bsc, n_differences=1):
        self._model_file = 'models/feasible_controller_placement_with_delay.mod'
        super(FeasibleControllerPlacementWithDelay, self).__init__(self._model_file)
        self.network = network
        self.M = M
        self.bcc = bcc
        self.bsc = bsc
        self.n_differences = n_differences

    def load_data(self, placements = []):
        super(FeasibleControllerPlacementWithDelay, self).load_data()

        T = {i: s for i, s in enumerate(placements)}

        vertex_list = set(self.network.nodes)
        self.ampl.set['VERTICES'] = set(vertex_list)
        self.ampl.set['L'] = set_range(placements)
        self.ampl.set['T'] = T
        self.ampl.set['U'] = pairs_exceeding_bcc_delay(self.network, self.bcc)
        self.ampl.set['W'] = nodes_within_bsc_delay(self.network, self.bsc)

        # print(f'M = {self.M}, m = {self.n_differences}')
        self.ampl.param['big_M'] = self.M
        self.ampl.param['small_m'] = self.n_differences

    def report(self):
        if not self._solved:
            self.solve()

        var_s = self.ampl.get_variables()['s']
        s_star = values_to_list(var_s)
        return {
            's': s_star,
        }

class MostDangerousKNodeAttack(MathematicalModel):
    def __init__(self, network, K):
        self._model_file = 'models/most_dangerous_k_node_attack.mod'
        super().__init__(self._model_file)
        self.K = K
        self.network = network

    def report(self):
        if not self._solved:
            self.solve()

        obj = self.ampl.getObjective('Payoff').value()
        return obj

    def generate_k(self, k):
        attacks = []

        for i in range(k):
            self.load_data(attacks)
            self.solve()
            # Get the current generated attack
            a = self.get_result_array('a')
            obj = self.ampl.getObjective('Payoff').value()
            assert(sum(a) == self.K)
            attacks.append(one_indices(a))
            print(f'objective: {obj}, attack = {one_indices(a)}')
        return attacks

    def load_data(self, attacks = []):
        super().load_data()

        nv = list(self.network.nodes)
        edges = self.network.edges

        def is_edge(x):
            x1, x2 = x
            return (x1, x2) in edges or (x2, x1) in edges

        adj_dict = {
            v: list(adjdict.keys()) for v, adjdict in self.network.adjacency()
        }

        h = list(filter(lambda x: x[0] < x[1], product(nv, nv)))
        he = list(filter(is_edge, h))

        deg = {
            v: self.network.degree(v) for v in nv
        }

        self.ampl.set['VERTICES'] = nv
        self.ampl.set['ATTACK_IDX'] = set_range(attacks)
        self.ampl.set['ATTACKS'] = to_set_ids(attacks)
        self.ampl.set['ADJ'] = adj_dict
        self.ampl.set['H'] = h
        self.ampl.set['H_EDGES'] = he

        self.ampl.param['deg'] = deg
        self.ampl.param['K'] = self.K


class FeasibleControllerPlacementWithDelayAndBC(MathematicalModel):
    def __init__(self, network, P, B, bcc, bsc, n_differences=1):
        self._model_file = 'models/feasible_controller_placement_with_delay_and_bc.mod'
        super().__init__(self._model_file)
        self.network = network
        self.P = P
        self.B  =B
        self.bcc = bcc
        self.bsc = bsc
        self.n_differences = n_differences

    def load_data(self, primary_controllers = [], backup_controllers = []):
        assert(len(primary_controllers) == len(backup_controllers))

        super().load_data()

        T_primary = {i: s for i, s in enumerate(primary_controllers)}
        T_backup = {i: s for i, s in enumerate(backup_controllers)}

        vertex_list = set(self.network.nodes)
        self.ampl.set['VERTICES'] = set(vertex_list)
        self.ampl.set['L'] = set_range(primary_controllers)
        self.ampl.set['T_primary'] = T_primary
        self.ampl.set['T_backup'] = T_backup
        self.ampl.set['U'] = pairs_exceeding_bcc_delay(self.network, self.bcc)
        self.ampl.set['W'] = nodes_within_bsc_delay(self.network, self.bsc)

        # print(f'M = {self.M}, m = {self.n_differences}')
        self.ampl.param['m'] = self.n_differences

        # Get the lowest
        self.ampl.param['P'] = self.P[0] if isinstance(self.P, tuple) else self.P
        self.ampl.param['B'] = self.B[0] if isinstance(self.B, tuple) else self.B


    def report(self):
        if not self._solved:
            self.solve()

        var_x = self.ampl.get_variables()['x']
        x_star = values_to_list(var_x)
        var_y = self.ampl.get_variables()['y']
        y_star = values_to_list(var_y)
        return {
            'x': x_star,
            'y': y_star,
        }

if __name__ == '__main__':
    from network import Network
    from algorithm import one_indices
    n = Network('cost266')

    p = AttackGenerationPricingProblem(n, (3, 6), budget=13, costs='degree')
    p.load_data([[4, 5, 24, 27, 30]], [1.0])
    r = p.report()
    a_star = one_indices(r['a*'])
    
    print(a_star)
    # primaries = []
    # backups = []

    # p = FeasibleControllerPlacementWithDelayAndBC(n, 3, 2, 3, 2)
    # p.load_data()
    # r = p.report()
    # prim = one_indices(r['y'])
    # backup = one_indices(r['x'])
    # print(prim, backup)

    # primaries.append(prim)
    # backups.append(backup)

    # for _ in range(3):
    #     p.load_data(primaries, backups)
    #     r = p.report()
    #     prim = one_indices(r['y'])
    #     backup = one_indices(r['x'])
    #     print(prim, backup)
    #     primaries.append(prim)
    #     backups.append(backup)
    
    # p = ControllerPlacementPricingProblemWithDelayAndBC(n, 1, 2, 3, 2)
    # p.load_data([[4, 5], [2, 8]], [0.4, 0.6])
    # print(p.report())
