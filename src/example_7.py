"""
Author: Gao Shuhua
Date: 2019/10/28
"""

from algorithm.bool_net import BooleanNetwork
from algorithm.fixed_destination_optimal_control import TimeVariantCostSolver


if __name__ == '__main__':
    L = [8, 4, 7, 7, 6, 2, 5, 5, 7, 3, 8, 8, 5, 1, 6, 6, 4, 8, 7, 7, 2, 6, 5, 5, 3, 7, 8, 8, 1, 5, 6, 6]
    m = 2
    n = 3
    T = 4
    x0 = 1
    Omega = {6}
    Cx = {1, 2, 3, 4, 5, 6, 7}


    def Cu(x):
        if x == 6:
            return {3, 4}
        else:
            return {1, 3, 4}


    def g(x, u, t):
        cu = [None, 2, 3, 1, 5 + t]
        return cu[u]


    def h(x, t):
        cx = [None, 3, 2 * t, 4, 0, 1, 5 + t, 6, 0]
        return cx[x]


    net = BooleanNetwork(L, m, n)
    solver = TimeVariantCostSolver(net, x0, Cx, Cu, Omega, h, g)
    J_star, u_star, s_star = solver.solve()
    print(f'J* = {J_star}')
    print(f'u* = {u_star}')
    print(f's* = {s_star}')

