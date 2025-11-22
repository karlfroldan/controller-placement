from amplpy import AMPL

def values_to_list(vs):
    return [x for _, x in vs.to_list()]

def to_set_ids(lst):
    """Convert the list into a dictionary of ID to item pair"""
    return { idx: x for idx, x in enumerate(lst) }

class MathematicalModel:
    def __init__(self, model_file):

        assert(model_file)
        self.ampl = AMPL()
        self._solved = False
        self._model_file = model_file
        self._model_loaded = False

    def _load_model(self, model_file = None):
        if not self._model_loaded:
            f = model_file if model_file else self._model_file
            self.ampl.read(f)
            self._model_loaded = True
    
    def solve(self, solver='cplex', model_file = None):
        self._load_model(model_file)
        # self.ampl.read(self._model_file)
        self.load_data()
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
    def __init__(self, network, placements, attacks):
        self._model_file = 'models/master_problem_mixed_strategies_main.mod';

        super(MixedStrategyMasterProblem, self).__init__(self._model_file)
        self._name = 'MixedStrategyMasterProblem'

        self.network = network
        self.placements = placements
        self.attacks = attacks
        

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
        v_star = self.ampl.get_variable('y').value()
        return {
            'V*': v_star,
            'q': qs,
            'p': ps,
        }

    def load_data(self):
        self._load_model()
        # Convert them to IDs so we can pass them as sets.
        placement_mapping = to_set_ids(self.placements)
        attack_mapping = to_set_ids(self.attacks)

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
    def __init__(self, network, placements, q_star, K, eps = 1e-9):
        self._model_file = 'models/attack_generation_pricing_problem.mod'
        super(AttackGenerationPricingProblem, self).__init__(self._model_file)
        self.network = network

        self._name = 'AttackGenerationPricingProblem'

        new_placements = []
        new_q_star = []

        # Do not consider placements that have 0 probability.
        for p, q in zip(placements, q_star):
            if q >= eps:
                new_placements.append(p)
                new_q_star.append(q)
        self.placements = new_placements
        self.q_star = new_q_star
        self.K = K

    def report(self):
        if not self._solved:
            self.solve()

        var_a = self.ampl.get_variables()['a']
        a_star = values_to_list(var_a)
        return {
            'a*': a_star,
        }
        

    def load_data(self):
        self._load_model()
        # Convert them to IDs
        s_prime = self.placements
        # s_prime = [i for i in range(len(self.placements))]
        placement_mapping = to_set_ids(s_prime)
        vertex_set = list(self.network.nodes)
        # vertex_set = list(self.network.nodes())
        edge_set = list(self.network.edges)

        v_s = placement_mapping

        self.ampl.set['V'] = vertex_set
        self.ampl.set['S_PRIME'] = list(placement_mapping.keys())
        self.ampl.set['E'] = edge_set
        self.ampl.set['V_S'] = v_s

        self.ampl.param['K'] = self.K
        self.ampl.param['q_star'] = self.q_star

class ControllerPlacementPricingProblem(MathematicalModel):
    def __init__(self, network, attacks, p_star, M, eps=1e-9):
        self._model_file = 'models/controller_placement_pricing_problem.mod'
        super(ControllerPlacementPricingProblem, self).__init__(self._model_file)

        self._name = 'ControllerPlacementPricingProblem'
        self.network = network
        self.attacks = attacks
        self.p_star = p_star
        self.M = M
        self.eps = eps

    def report(self):
        if not self._solved:
            self.solve()
        var_s = self.ampl.get_variables()['s']
        s_star = values_to_list(var_s)
        # Y_star = values_to_list(self.ampl.get_variables()['Y'])

        V_star = self.ampl.get_objectives()['OperatorPayoff'].value()
        print(V_star)
        return {
            's*': s_star,
            'V*': V_star,
        }

    def load_data(self):
        self._load_model()
        vertex_list = list(self.network.nodes)

        # Key => attack set pair
        attack_dict = to_set_ids(self.attacks)

        self.ampl.set['VERTICES'] = vertex_list
        self.ampl.set['ATTACKS'] = attack_dict.keys()

        
        V_A = {}
        for i, a in attack_dict.items():
            # vs = self.network.remaining_nodes(a)
            V_A[i] = a
        self.ampl.set['V_A'] = V_A

        component_ids = {}
        C_A = {}
        for a, attack in attack_dict.items():
            components = self.network.attack(attack)
            component_ids[a] = list(range(len(components)))
            for i, c in enumerate(components):
                C_A[(a, i)] = list(c)

        print(component_ids)
        self.ampl.set['COMPONENT_IDS'] = component_ids
        self.ampl.set['C_A'] = C_A

        self.ampl.param['M'] = self.M
        self.ampl.param['p_star'] = self.p_star
