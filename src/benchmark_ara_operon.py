"""
The ara operon network in E. coli.

Author: Gao Shuhua
Date: 2019/10/28
"""
from algorithm.bool_net import BooleanNetwork, STG
from algorithm.fixed_time_optimal_control import FixedTimeSolver
from algorithm.fixed_destination_optimal_control import TimeInvariantCostSolver
from algorithm.utils import read_network
from algorithm.existing_work import *
import time


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


def inverse_map(i: int, n: int):
    """
    Accumulative STP of logical variables is bijective.
    Given a result i (\delta_{2^n}^i), find the corresponding logical values.

    :return a list of 0/1
    """
    r = []
    while n > 0:
        if i % 2 == 0:
            r.append(0)
            i = i // 2
        else:
            r.append(1)
            i = (i + 1) // 2
        n = n - 1
    r.reverse()
    return r


def g(x, u, n, m):
    X = inverse_map(x, n)
    U = inverse_map(u, m)
    A = [0, 16, 40, 44, 28, 28, 28, 48, 44]
    B = [0, 48, 28, 48]
    return sum(a * x for a, x in zip(A, X)) + sum(b * u for b, u in zip(B, U))


def run_task1(n, m, L, x0, xd, T, which_solver, solver_name: str):
    print(f'- Running [{solver_name}] for Task 1')
    ts = time.time()
    solver = which_solver(L, m, n, lambda x, u: g(x, u, n, m))
    J_star = solver.solve(x0, xd, T)
    print(f'\tFinished. J* = {J_star}. Time elapsed (s): {time.time() -ts: .2f}')


def run_task2(n, m, L, x0, xd, which_solver, solver_name: str):
    print(f'- Running [{solver_name}] for Task 2')
    ts = time.time()
    solver = which_solver(L, m, n, lambda x, u: 1)
    T_star = solver.solve(x0, xd)
    print(f'\tFinished. T* = {T_star}. Time elapsed (s): {time.time() -ts: .2f}')


def visualize(x0, xd, T):
    n, m, L = read_network('ara_operon.txt')
    data = {}
    # STG
    net = BooleanNetwork(L, m, n)
    stg = STG(net, x0, None, None)
    for i in stg.vertices:
        for j in stg.adj_list[i]:
            data[(i, j)] = '->'
    # min energy
    solver = FixedTimeSolver(net, x0, T, None, None, {xd}, None, lambda x, u, t: g(x, u, n, m))
    J_star, u_star, s_star = solver.solve()
    for t in range(len(s_star) - 1):
        data[(s_star[t], s_star[t + 1])] = 'me'
    # min time
    solver = TimeInvariantCostSolver(net, x0, None, None, {xd}, None, lambda x, u, t: 1)
    J_star, u_star, s_star = solver.solve()
    for t in range(len(s_star) - 1):
        if data[(s_star[t], s_star[t + 1])] == 'me':
            data[(s_star[t], s_star[t + 1])] = 'mm'
        else:
            data[(s_star[t], s_star[t + 1])] = 'mt'
    with open('benchmark_ara_operon.sif', 'w') as f:
        for (i, j), sym in data.items():
            f.write(f'{i} {sym} {j}\n')


if __name__ == '__main__':
    n, m, L = read_network('ara_operon.txt')
    x0 = 9
    xd = 410
    print('Task 1: minimum-energy control')
    T = 10
    # visualize(x0, xd, T)
    run_task1(n, m, L, x0, xd, T, ProposedSolverForTask1, 'Proposed approach')
    run_task1(n, m, L, x0, xd, T, EttoreFornasiniSolver, 'Ettore Fornasini etc.')
    run_task1(n, m, L, x0, xd, T, CuiXingbangSolver, 'Cui Xingbang etc.')
    run_task1(n, m, L, x0, xd, T, ZhuQunxiSolver, 'Zhu Qunxi etc.')
    run_task1(n, m, L, x0, xd, T, LiFangfeiSolver, 'Li Fangfei etc.')

    print('Task 2: minimum-time control')
    run_task2(n, m, L, x0, xd, ProposedSolverForTask2, 'Proposed approach')
    run_task2(n, m, L, x0, xd, DmitriyLaschovSolver, 'Dmitriy Laschov etc.')
    run_task2(n, m, L, x0, xd, LiFangfeiSolver, 'Li Fangfei etc.')  #note: a very long time