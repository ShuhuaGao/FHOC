"""
The T-LGL network in leukemia

Author: Gao Shuhua
Date: 2020/3/19
"""
from algorithm.bool_net import BooleanNetwork, STG
from algorithm.fixed_time_optimal_control import FixedTimeSolver
from algorithm.fixed_destination_optimal_control import TimeInvariantCostSolver
from algorithm.utils import read_network
from algorithm.existing_work import *
import time
import numpy as np


np.random.seed(0)   # for reproduction
# the weight matrix in specifying the stage cost
N = 65536
M = 8
Q = np.random.uniform(0, 5, N)
R = np.random.uniform(0, 10, M)


class ProposedSolverForTask1(AbstractSolver):
    """
    A simple wrapper of our graph-theoretical approach.
    """
    def __init__(self, L: List[int], m: int, n: int, g: Callable[[int, int], float]):
        super().__init__(L, m, n, g)
        self.net = BooleanNetwork(L, m, n)

    def solve(self, x0: int, xd: int, T: int=None) -> float:
        if T is not None:
            solver = FixedTimeSolver(self.net, x0, T, None, None, {xd}, None, lambda x, u, t: self.g(x, u))
            J_star, u_star, s_star = solver.solve()
        else: # time is not specified
            net = BooleanNetwork(L, m, n)
            solver = TimeInvariantCostSolver(net, x0, None, None, {xd}, None, lambda x, u, t: self.g(x, u))
            J_star, u_star, s_star = solver.solve()
        J = 0
        for ut, xt in zip(u_star, s_star[:T]):
            J = J + self.g(xt, ut)
        assert J == J_star
        return J_star


class ProposedSolverForTask2(AbstractSolver):
    """
        A simple wrapper of our graph-theoretical approach.
    """
    def __init__(self, L: List[int], m: int, n: int, g: Callable[[int, int], float]):
        super().__init__(L, m, n, g)
        self.net = BooleanNetwork(L, m, n)

    def solve(self, x0: int, xd: int, T: int=None) -> int:
        solver = TimeInvariantCostSolver(self.net, x0, None, None, {xd}, None, lambda x, u, t: 1)
        J_star, u_star, s_star = solver.solve()
        return J_star


def g(x: int, u: int) -> float:
    """
    Stage cost function

    :param x: i in a state $\delta_N^i$
    :param u: k in a control $\delta_M^k$
    :return: cost
    """
    return Q[x - 1] + R[u - 1]


def run_task1(n, m, L, x0, xd, T, which_solver, solver_name: str):
    print(f'- Running [{solver_name}] for Task 1')
    ts = time.time()
    solver = which_solver(L, m, n, g)
    J_star = solver.solve(x0, xd, T)
    print(f'\tFinished. J* = {J_star}. Time elapsed (s): {time.time() -ts: .2f}')


def run_task2(n, m, L, x0, xd, which_solver, solver_name: str):
    print(f'- Running [{solver_name}] for Task 2')
    ts = time.time()
    solver = which_solver(L, m, n, lambda x, u: 1)
    T_star = solver.solve(x0, xd)
    print(f'\tFinished. T* = {T_star}. Time elapsed (s): {time.time() -ts: .2f}')


if __name__ == '__main__':
    n, m, L = read_network('T_LGL.txt')
    x0 = 58834
    xd = 65535
    print('Task 1: minimum-energy control')
    T = 10

    run_task1(n, m, L, x0, xd, T, ProposedSolverForTask1, 'Proposed approach')
    run_task1(n, m, L, x0, xd, T, CuiXingbangSolver, 'Cui Xingbang etc.')
    run_task1(n, m, L, x0, xd, T, ZhuQunxiSolver, 'Zhu Qunxi etc.')
    run_task1(n, m, L, x0, xd, T, EttoreFornasiniSolver, 'Ettore Fornasini etc.')
    run_task1(n, m, L, x0, xd, T, LiFangfeiSolver, 'Li Fangfei etc.')



    print('Task 2: minimum-time control')
    run_task2(n, m, L, x0, xd, ProposedSolverForTask2, 'Proposed approach')
    run_task2(n, m, L, x0, xd, LiFangfeiSolver, 'Li Fangfei etc.')  # note: a very long time
    run_task2(n, m, L, x0, xd, DmitriyLaschovSolver, 'Dmitriy Laschov etc.')


