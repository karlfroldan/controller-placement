import time
import numpy as np

from network import Network

from model import *
from utils import make_ones, one_indices

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

class InitialGenerationWithDelay:
    def __init__(self, network, M, K, bsc, bcc, n_generations):
        self.M = M
        self.K = K
        self.network = network
        self.placement_gen = FeasibleControllerPlacementWithDelay(network, M, bcc, bsc)
        self.n_placement_gens = n_generations

    def initialize(self):
        V = len(list(self.network.nodes))
        attack = one_indices(make_ones(V, self.K))
        placements = []
        for i in range(self.n_placement_gens):
            if i == 0:
                self.placement_gen.load_data()
            else:
                self.placement_gen.load_data(placements = placements)


            self.placement_gen.solve()
            r = self.placement_gen.report()
            if one_indices(r['s']) == []:
                break
            s = one_indices(r['s'])
            placements.append(s)

        return {
            'placement': placements,
            'attack': [attack],
        }

class InitialGenerationWithDelayAndBC:
    def __init__(self, network, C, P, B, K, bsc, bcc, n_generations):
        self.P = P
        self.B = B
        self.K = K
        self.network = network

        self.placement_gen = FeasibleControllerPlacementWithDelayAndBC(network, P, B, C, bcc, bsc)
        self.n_placement_gens = n_generations

    def initialize(self):
        V = len(list(self.network.nodes))
        K = self.K[0] if isinstance(self.K, tuple) else self.K
        attack = one_indices(make_ones(V, K))

        primaries = []
        backups = []

        for i in range(self.n_placement_gens):
            # if i == 0:
            #     self.placement_gen.load_data()
            # else:
            #     self.placement_gen.load_data(primaries, backups)
            self.placement_gen.load_data(primaries, backups)
            self.placement_gen.solve()
            r = self.placement_gen.report()
            if one_indices(r['y']) == [] and one_indices(r['x']) == []:
                break

            prim = one_indices(r['y'])
            backup = one_indices(r['x'])

            primaries.append(prim)
            backups.append(backup)

        print(primaries)
        print(backups)

        return {
            'primary_controllers': primaries,
            'backup_controllers': backups,
            'attack': [attack],
        }

class ControllerPlacementOptimization:
    def __init__(self, network, M, K, eps = 1e-9):
        self.M = M
        self.K = K
        self.network = network
        self.V = len(list(network.nodes))

        self.initializer_model = InitialGeneration(network, M, K)

        self.attack_gen_model = AttackGenerationPricingProblem(network, K)
        self.controller_gen_model = ControllerPlacementPricingProblem(network, M)

        self.master_model = MixedStrategyMasterProblem(network)

        self.eps = eps

    def run(self):
        initial = self.initializer_model.initialize()
        placements = initial['placement']
        attacks = initial['attack']

        print(f'Initial placements: {placements}')
        print(f'Initial attacks: {attacks}')
        
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
        # mp = MixedStrategyMasterProblem(self.network, placements, attacks)
        self.master_model.load_data(placements, attacks)
        self.master_model.solve()
        mp_values = self.master_model.report()

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
        self.controller_gen_model.load_data(attacks, p_star)
        self.controller_gen_model.solve()
        return self.controller_gen_model.report()
        
    def _make_attack_generation_pricing_problem(self, placements, q_star):
        self.attack_gen_model.load_data(placements, q_star)
        self.attack_gen_model.solve()
        
        return self.attack_gen_model.report()

    def _one_round(self, placements, attacks, p_star, q_star, x_star, y_star):
        updated_placement = False
        # Solve the placement generation problem CP[ATTACKS, p*] to get placement
        # s'.
        pp = self._make_controller_pricing_problem(attacks, p_star)
        s_star = one_indices(pp['s*'])
        
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
        a_star = one_indices(agp['a*'])

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
    def __init__(self, network, M, K, bsc, bcc, initialized_placements = 1):
        super().__init__(network, M, K)
        self.controller_gen_model = ControllerPlacementPricingProblemWithDelay(network, M, bsc, bcc)
        self.initializer_model = InitialGenerationWithDelay(network, M, K, bsc, bcc, initialized_placements)
        self.bsc = bsc
        self.bcc = bcc

class ControllerOptimizationWithDelayAndBC(ControllerPlacementOptimization):
    def __init__(self, network, C, P, B, K, bsc, bcc, initialized_placements = 1, attacker_budget=0, attacker_costs=None):
        super().__init__(network, 0, K)
        
        self.P = P
        self.B = B
        self.K = K

        self.controller_gen_model = ControllerPlacementPricingProblemWithDelayAndBC(network, C, P, B, bsc, bcc)
        self.attack_gen_model = AttackGenerationPricingProblem(network, K, budget=attacker_budget, costs=attacker_costs)
        self.initializer_model = InitialGenerationWithDelayAndBC(network, C, P, B, K, bsc, bcc, initialized_placements)
        self.bsc = bsc
        self.bcc = bcc

    def run(self):
        initial = self.initializer_model.initialize()
        primary_placements = initial['primary_controllers']
        backup_placements = initial['backup_controllers']
        attacks = initial['attack']

        print(f'Initial placements: {primary_placements}, {backup_placements}')
        print(f'Initial attacks: {attacks}')
        
        # Solve the master problem first.
        mp_out = self._solve_master_problem((primary_placements, backup_placements), attacks)

        placements = (primary_placements, backup_placements)

        total_placement_time = 0.0
        total_attack_time = 0.0

        while True:
            mp_out = self._one_round(placements, attacks, mp_out['p*'], mp_out['q*'], mp_out['x*'], mp_out['y*'])
            total_placement_time += mp_out['placement_time']
            total_attack_time += mp_out['attack_time']
            
            if not mp_out['updated']:
                break;

        return {
            'V*': mp_out['x*'],
            'p*': mp_out['p*'],
            'q*': mp_out['q*'],
            'primary_placements': mp_out['primary_placements'],
            'backup_placements': mp_out['backup_placements'],
            'attacks': mp_out['attacks'],
            'placement_time': total_placement_time,
            'attack_time': total_attack_time,
        }        

    def _one_round(self, placements, attacks, p_star, q_star, x_star, y_star):
        updated_placement = False
        # Solve the placement generation problem CP[ATTACKS, p*] to get placement
        # s'.
        t_begin = time.time()
        pp = self._make_controller_pricing_problem(attacks, p_star)
        t_end = time.time()
        pc_time = t_end - t_begin

        primary_placement = one_indices(pp['y*'])
        backup_placement = one_indices(pp['x*'])

        
        lhs = sum([
            len(self.network.surviving_nodes(primary_placement, a, backup_controllers=backup_placement)) * p_a
            for a, p_a in zip(attacks, p_star)
        ])

        rhs = x_star

        if lhs > rhs + self.eps:
            placements[0].append(primary_placement)
            placements[1].append(backup_placement)
            updated_placement = True

        # Solve the master problem again
        mp_out = self._solve_master_problem(placements, attacks)
        q_star = mp_out['q*']

        # Solve the attack generation problem
        # For the attack generation problem, we might only want to care about the primary controllers
        t_begin = time.time()
        agp = self._make_attack_generation_pricing_problem(placements, q_star)
        t_end = time.time()
        ag_time = t_end - t_begin
        a_star = one_indices(agp['a*'])

        rhs = sum([
            len(self.network.surviving_nodes(s1, a_star, backup_controllers=s2)) * q_a
            for s1, s2, q_a in zip(placements[0], placements[1], q_star)
        ])

        lhs = y_star
        if lhs > rhs + self.eps:
            attacks.append(a_star)
            updated_placement = True

        mp_out = self._solve_master_problem(placements, attacks)

        return {
            'updated': updated_placement,
            'primary_placements': placements[0],
            'backup_placements': placements[1],
            'attacks': attacks,
            'x*': mp_out['x*'],
            'y*': mp_out['y*'],
            'q*': mp_out['q*'],
            'p*': mp_out['p*'],
            'placement_time': pc_time,
            'attack_time': ag_time,
        }

class PureControllerPlacementGeneration:
    def __init__(self, network, M, K):
        self.M = M
        self.K = K
        self.naop = NAOP(network)
        self.cpop = CPOP(network)
        self.V = len(list(network.nodes))

    def run(self, attacks = []):
        # Generate a random M-node controller placement s*
        s_star = one_indices(make_ones(self.V, self.M))
        Y_star = self.V
        Z_star = 0

        n_iterations = 0

        total_cpop_time = 0
        total_naop_time = 0

        while True:
            self.naop.load_data([s_star], self.K)

            begin_time = time.time()
            self.naop.solve()
            end_time = time.time()

            total_naop_time += end_time - begin_time
            
            r = self.naop.report()
            a_star = one_indices(r['a'])
            Z_star = r['Z']

            if Z_star >= Y_star:
                break

            attacks.append(a_star)

            # Solve CPOP to get best placement s*
            self.cpop.load_data(attacks, self.M)

            begin_time = time.time()
            self.cpop.solve()
            end_time = time.time()

            total_cpop_time = end_time - begin_time

            r = self.cpop.report()
            s_star = one_indices(r['s'])
            Y_star = r['Y']

            # Y* is equal to
            # V(s*, a(s*)) = max_{s in S(M)} min(a in A) V(s, a)
            n_iterations += 1
            
        print(f'Done in {n_iterations} iterations')
        print(f'Surviving Nodes: {Y_star}')
        return {
            'attacks': attacks,
            's*': s_star,
            'Y*': Y_star,

            'total_cpop_time': total_cpop_time,
            'total_naop_time': total_naop_time,
        }

class PureAttackGeneration:
    def __init__(self, network, M, K):
        self.M = M
        self.K = K
        self.naop = NAOP(network)
        self.cpop = CPOP(network)
        self.V = len(list(network.nodes))

    def run(self, placements = []):
        a_star = one_indices(make_ones(self.V, self.K))
        Z_star = 0
        Y_star = self.V

        n_iterations = 0

        while True:
            self.cpop.load_data([a_star], self.M)
            self.cpop.solve()
            r = self.cpop.report()
            s_star = one_indices(r['s'])
            Y_star = r['Y']

            if Y_star <= Z_star:
                break
            placements.append(s_star)
            
            self.naop.load_data(placements, self.K)
            self.naop.solve()
            r = self.naop.report()
            a_star = one_indices(r['a'])
            Z_star = r['Z']
            n_iterations += 1
            
        print(f'Done in {n_iterations} iterations')
        print(f'Surviving Nodes: {Z_star}')
        return {
            'placements': placements,
            'a*': a_star,
            'Z*': Z_star,
        }

if __name__ == '__main__':
    n = Network('cost266')
    eps = 1e-9
    M = 6
    K = 3

    # p = ControllerOptimizationWithDelay(n, M, K, 1529.28, 1500)
    p = ControllerPlacementOptimization(n, M, K)
    r = p.run()
    print(r)
