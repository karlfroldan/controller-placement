module OptimUtils

using JuMP
import Iterators
import CPLEX

using Graphs

# function directed_node_pairs(g :: SimpleGraph)
#     V = nv(g)

#     V² = filter(v -> v[1] < v[2], Iterators.product(1:V, 1:V))

#     H = Set{Tuple{Int, Int}}()
#     H_prime = Set{Tuple{Int, Int}}()

#     for (v, w) in V²
#         for e in edges(g)
#             if src(e) == v && dst(e) == w && v < w
#                 push!(H, 
#             end
#         end
#     end
# end

# function generate_topological_attack(g :: SimpleGraph, K :: Int)
#     m = Model(CPLEX.Optimizer)
# end

end
