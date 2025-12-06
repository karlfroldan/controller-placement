from model import AttackGenerationPricingProblem, ControllerPlacementPricingProblem
from network import Network
from algorithm import one_indices, ControllerPlacementOptimization, PureControllerPlacementGeneration, PureAttackGeneration, FeasibleControllerPlacementWithDelay, ControllerOptimizationWithDelay

def add1(x):
    return x + 1

# p = PureControllerPlacementGeneration(n, 6, 4)
n = Network('cost266')
p = ControllerOptimizationWithDelay(n, 6, 4, 2000, 2000)
r = p.run()
# print(f'Y*: {r["Y*"]}, Attack length: {len(r["attacks"])}')
# cpop_time = r['total_cpop_time']
# naop_time = r['total_naop_time']
# print(f'cpop time: {cpop_time}, naop_time: {naop_time}')
print(f'V*: {r["V*"]}, Attack length: {len(r["attacks"])}')

print('-=-=-=-=-= PLACEMENTS -=-=-=-=-=')
for q, s in zip(r['q*'], r['placements']):
    print(f'{q:.2f} => {s}')

print('-=-=-=-=-= ATTACKS -=-=-=-=-=')
for p, a in zip(r['p*'], r['attacks']):
    print(f'{p:.2f} => {a}')

# p = FeasibleControllerPlacementWithDelay(n, 6, 2000, 2000)
# p.load_data()
# p.solve()
# r = p.report()
# print(r)
