"""
Author: Gao Shuhua
Date: 2019/10/27
"""

from algorithm.bool_net import BooleanNetwork
from algorithm.fixed_destination_optimal_control import STGPlus, STG, TimeInvariantCostSolver

if __name__ == '__main__':
    L = [8, 4, 7, 7, 6, 2, 5, 5, 7, 3, 8, 8, 5, 1, 6, 6, 4, 8, 7, 7, 2, 6, 5, 5, 3, 7, 8, 8, 1, 5, 6, 6]
    m = 2
    n = 3
    T = 4
    x0 = 7
    Omega = {3, 4}
    Cx = {1, 2, 3, 4, 5, 6, 7}


    def Cu(x):
        if x == 6:
            return {3, 4}
        else:
            return {1, 3, 4}


    def g(x, u, t):
        cx = [None, 2, 5, 1, 4, 1, 3, 6, 0]
        cu = [None, 0, 3, 1, 4]
        return cx[x] + cu[u]


    def h(x, K):
        return 0


    net = BooleanNetwork(L, m, n)


    solver = TimeInvariantCostSolver(net, x0, Cx, Cu, Omega, h, g)
    J_star, u_star, s_star = solver.solve()
    print(f'J^* = {J_star}')
    print(f'u^* = {u_star}')
    print(f's^* = {s_star}')