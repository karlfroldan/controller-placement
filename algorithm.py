import numpy as np

from dataloader import Network

from models import *

def make_ones(m, n):
    """Create a list of size M with N ones"""
    assert(m >= n)
    ones = np.ones(n, dtype='int')
    zeros = np.zeros(m - n, dtype='int')

    rng = np.random.default_rng()
    arr = np.concatenate((ones, zeros))
    rng.shuffle(arr)
    return arr

def one_indices(xs):
    # add1 = lambda x: x + 1

    if isinstance(xs, list):
        xs = np.array(xs)
    return np.where(xs == 1)[0].tolist()

class InitialGeneration:
    """
    Initial generation of the attack and placement list class.
    This class simply creates a basic placement and an attack.
    """
    def __init__(self, network, M, K):
        self.M = M
        self.K = K
        self.network = network

    def initialize(self):
        V = len(list(self.network.nodes))
        return {
            'placement': [one_indices(make_ones(V, self.M))],
            'attack': [one_indices(make_ones(V, self.K))],
        }

class ControllerPlacementOptimization:
    def __init__(self, network, M, K, eps = 1e-9):
        self.M = M
        self.K = K
        self.network = network
        self.V = len(list(network.nodes))

        assert(self.V >= M)
        assert(self.V >= K)

        self.initializer_model = InitialGeneration(network, M, K)

        self.attack_gen_model = AttackGenerationPricingProblem
        self.controller_gen_model = ControllerPlacementPricingProblem

        self.eps = eps

    def run(self):
        initial = self.initializer_model.initialize()
        placements = initial['placement']
        attacks = initial['attack']
        
        # Solve the master problem first.
        mp_out = self._solve_master_problem(placements, attacks)

        while True:
            mp_out = self._one_round(placements, attacks, mp_out['p*'], mp_out['q*'], mp_out['x*'], mp_out['y*'])
            if not mp_out['updated']:
                break;

        return {
            'V*': mp_out['x*'],
            'p*': mp_out['p*'],
            'q*': mp_out['q*'],
            'placements': mp_out['placements'],
            'attacks': mp_out['attacks'],
        }


    def _solve_master_problem(self, placements, attacks):
        mp = MixedStrategyMasterProblem(self.network, placements, attacks)
        mp.solve()
        mp_values = mp.report()

        p_star = mp_values['p']
        q_star = mp_values['q']

        x_star = mp_values['x']
        y_star = mp_values['y']

        return {
            'p*': p_star,
            'q*': q_star,
            'x*': x_star,
            'y*': y_star,
        }

    def _make_controller_pricing_problem(self, attacks, p_star):
        return self.controller_gen_model(self.network, attacks, p_star, self.M)
        
    def _make_attack_generation_pricing_problem(self, placements, q_star):
        return self.attack_gen_model(self.network, placements, q_star, self.K)

    def _one_round(self, placements, attacks, p_star, q_star, x_star, y_star):
        updated_placement = False
        # Solve the placement generation problem CP[ATTACKS, p*] to get placement
        # s'.
        pp = self._make_controller_pricing_problem(attacks, p_star)
        pp_values = pp.report()
        s_star = one_indices(pp_values['s*'])

        lhs = sum([
            len(self.network.surviving_nodes(s_star, a)) * p_a
            for a, p_a in zip(attacks, p_star)
        ])

        rhs = x_star

        # print(lhs, rhs)
        if lhs > rhs + self.eps:
            placements.append(s_star)
            updated_placement = True

        # Solve the master problem again
        mp_out = self._solve_master_problem(placements, attacks)
        q_star = mp_out['q*']

        # Solve the attack generation problem
        agp = self._make_attack_generation_pricing_problem(placements, q_star)
        agp_values = agp.report()
        a_star = one_indices(agp_values['a*'])

        rhs = sum([
            len(self.network.surviving_nodes(s, a_star)) * q_a
            for s, q_a in zip(placements, q_star)
        ])

        lhs = y_star
        if lhs > rhs + self.eps:
            attacks.append(a_star)
            updated_placement = True

        mp_out = self._solve_master_problem(placements, attacks)

        return {
            'updated': updated_placement,
            'placements': placements,
            'attacks': attacks,
            'x*': mp_out['x*'],
            'y*': mp_out['y*'],
            'q*': mp_out['q*'],
            'p*': mp_out['p*'],
        }

class ControllerOptimizationWithDelay(ControllerPlacementOptimization):
    def __init__(self, network, M, K, bsc, bcc):
        super(ControllerOptimizationWithDelay, self).__init__(network, M, K)
        self.controller_gen_model = ControllerPlacementPricingProblemWithDelay
        self.bsc = bsc
        self.bcc = bcc

    def _make_controller_pricing_problem(self, attacks, p_star):
        return self.controller_gen_model(self.network, attacks, p_star, self.M, self.bsc, self.bcc)

if __name__ == '__main__':
    n = Network('dognet')
    eps = 1e-9
    opt = ControllerOptimizationWithDelay(n, 2, 2, 3, 3)
    vals = opt.run()
    print(f'Placements')
    for s, q in zip(vals['placements'], vals['q*']):
        if q >= eps:
            print(f'{s} - {q:.2f}')
    print('Attacks')
    for a, p in zip(vals['attacks'], vals['p*']):
        if p >= eps:
            print(f'{a} - {p:.2f}')
    # print(opt.run())
    
