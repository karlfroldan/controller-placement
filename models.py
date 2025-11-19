from amplpy import AMPL

class MixedStrategyMasterProblem(AMPL):
    def __init__(self, network, placements, attacks):
        super(MixedStrategyMasterProblem, self).__init__()
        
        self.network = network
        self.placements = placements
        self.attacks = attacks

        self.solved = False

    def solve(self, solver='cplex'):
        self.read('models/master_problem_mixed_strategies_main.mod')
        self.load_data()
        self.option['solver'] = solver
        super(MixedStrategyMasterProblem, self).solve()
        self.solved = True

    def report(self):
        if not self.solved:
            self.solve()

        def values_to_list(vs):
            fst = lambda x: x[0]
            snd = lambda x: x[1]
            return list(map(snd, sorted(vs.to_list(), key=fst)))

        # This is the probability distribution of the attacker mixed strategy
        ps = self.get_constraint('BoundPayoff').get_values()
        ps = values_to_list(ps)

        qs = self.get_variables()['q']
        # Return the controller placement mixed-strategy.
        qs = values_to_list(qs)

        # payoff
        v_star = self.get_variable('y').value()
        return {
            'V*': v_star,
            'q': qs,
            'p': ps,
        }

    def load_data(self):
        # Convert them to IDs so we can pass them as sets.
        placement_mapping = { idx: p for idx, p in enumerate(self.placements) }
        attack_mapping = { idx: a for idx, a in enumerate(self.attacks) }

        self.set['S'] = sorted(placement_mapping.keys())
        self.set['A'] = sorted(attack_mapping.keys())

        v_data = {}

        for i, p in placement_mapping.items():
            for j, a in attack_mapping.items():
                surviving_nodes = len(self.network.surviving_nodes(p, a))
                assert(isinstance(surviving_nodes, int))

                v_data[(i, j)] = surviving_nodes

        self.param['V'] = v_data
