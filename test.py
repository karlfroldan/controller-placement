from model import *
from network import Network
from algorithm import *
from utils import one_indices

def add1(xs):
    list(map(lambda x: x + 1, xs))

# p = PureControllerPlacementGeneration(n, 6, 4)

n = Network('cost266')

pc_range = (5, 9)
bc_range = (2, 4)
a_range = (3, 9)
p = ControllerOptimizationWithDelayAndBC(n, pc_range, bc_range, a_range, 2000, 2000, attacker_budget=20, attacker_costs='degree')
r = p.run()

# print(f'Y*: {r["Y*"]}, Attack length: {len(r["attacks"])}')
# cpop_time = r['total_cpop_time']
# naop_time = r['total_naop_time']
# print(f'cpop time: {cpop_time}, naop_time: {naop_time}')
# print(f'V*: {r["V*"]}, Attack length: {len(r["attacks"])}')

print(f'TIME ELAPSED ({r["placement_time"]}, {r["attack_time"]}')

print('-=-=-=-=-= PLACEMENTS -=-=-=-=-=')
for q, s1, s2 in zip(r['q*'], r['primary_placements'], r['backup_placements']):
    if q > 1e-9:
        print(f'{q:.3f} => {s1}, backup = {s2}')

print('-=-=-=-=-= ATTACKS -=-=-=-=-=')
for p, a in zip(r['p*'], r['attacks']):
    if p >  1e-9:
        print(f'{p:.3f} => {a}')

print(f'\n\nPAYOFF: {r["V*"]}')

# p = FeasibleControllerPlacementWithDelay(n, 20, 2000, 2000)
# p.load_data()
# p.solve()
# r = p.report()
# print(r)
