"""
Author: Gao Shuhua
Date: 2019/10/26
"""
from algorithm.bool_net import BooleanNetwork
from algorithm.fixed_time_optimal_control import FixedTimeSolver


if __name__ == '__main__':
    L = [ 8, 4, 7, 7, 6, 2, 5, 5, 7, 3, 8, 8, 5, 1, 6, 6, 4, 8, 7, 7, 2, 6, 5, 5, 3, 7, 8, 8, 1, 5, 6, 6]
    m = 2
    n = 3
    T = 4
    x0 = 1
    Omega = {2, 6}
    Cx = {1, 2, 3, 4, 5, 6, 7}

    def Cu(x):
        if x == 6:
            return {3, 4}
        else:
            return {1, 3, 4}

    def g(x, u, t):
        c = [None, 2, 3, 1, 0]
        return c[u] + t

    def h(x):
        c = [None, 3, 5, 4, 0, 1, 3, 6, 0]
        return c[x]

    net = BooleanNetwork(L, m, n)
    solver = FixedTimeSolver(net, x0, T, Cx, Cu, Omega, h, g)
    J_star, u_star, s_star = solver.solve()
    print(f'J_T^* = {J_star}')
    print(f'u^* = {u_star}')
    print(f's^* = {s_star}')