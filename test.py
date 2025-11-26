from models import AttackGenerationPricingProblem, ControllerPlacementPricingProblem
from network import Network
from algorithm import one_indices, ControllerPlacementOptimization

def add1(x):
    return x + 1

n = Network('cost266')

p = ControllerPlacementOptimization(n, 13, 6)
r = p.run()
print(f'V*: {r["V*"]}, Attack length: {len(r["attacks"])}')
