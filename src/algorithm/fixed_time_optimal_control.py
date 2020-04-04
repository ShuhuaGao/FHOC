"""
Fixed-time optimal control
Refer to Algorithm 2 in the paper.

Author: Gao Shuhua
Date: 2019/10/24
"""
from typing import Callable, List, Dict, Iterable, Tuple
from operator import itemgetter
from collections import defaultdict
from .bool_net import BooleanNetwork


class FixedTimeSolver:
    """
    Solver for fixed-time optimal control
    """
    def __init__(self, net: BooleanNetwork, x0: int, T: int,
                 Cx: Iterable[int], Cu: Callable[[int], Iterable[int]], Omega: Iterable[int],
                 h: Callable[[int], float], g: Callable[[int, int, int], float]):
        """
        Initialization.

        :param net: the network
        :param x0: initial state
        :param T: horizon length
        :param Cx: state constraints, only states in Cx are allowed
        :param Cu: control constraints, Cu(i) gives the control that are allowed at i
        :param Omega: terminal constraint
        :param h: terminal cost
        :param g: stage cost
        """
        self.net = net
        self.Cx = Cx
        self.Cu = Cu
        self.x0 = x0
        self.T = T
        if isinstance(Omega, set):
            self.Omega = Omega
        else:
            self.Omega = set(Omega)
        if g is not None:
            self.g = g
        else:
            self.g = lambda x, u, t: 0
        if h is not None:
            self.h = h
        else:
            self.h = lambda x: 0

    def _build_TE_STG(self):
        """
        Get the TE-STG following a BFS procedure in a layer-by-layer manner.
        """
        _successors = {}
        predecessors = defaultdict(list)
        layers = [set() for _ in range(self.T + 2)] # layers of vertices
        layers[0].add(self.x0)
        for t in range(self.T - 1):
            for i in layers[t]:
                if i not in _successors:
                    _successors[i] = self.net.successors(i, self.Cx, self.Cu)
                layers[t + 1] = layers[t + 1] | _successors[i].keys()
                for j in _successors[i]:
                    ks = _successors[i][j]  # control inputs that attain transition from i to j
                    opt = min(((k, self.g(i, k, t)) for k in ks), key=itemgetter(1))
                    predecessors[(t+1, j)].append((i, opt[0], opt[1])) # (predecessor, control, weight)
        # the T layer, we need to consider the terminal constraint
        for i in layers[self.T - 1]:
            if i not in _successors:
                _successors[i] = self.net.successors(i, self.Cx, self.Cu)
            for j in (self.Omega & _successors[i].keys()):
                layers[self.T].add(j)
                ks = _successors[i][j]  # control inputs that attain transition from i to j
                opt = min(((k, self.g(i, k, self.T - 1)) for k in ks), key=itemgetter(1))
                predecessors[(self.T, j)].append((i, opt[0], opt[1]))
        # the pseudo-state layer
        layers[self.T + 1].add(0)
        for i in layers[self.T]:
            w = self.h(i)
            predecessors[(self.T + 1, 0)].append((i, None, w))
        # assert existence of solutions: Proposition 1
        assert layers[self.T], "Problem is not feasible! Check Proposition 1"
        # predecessors[(t, j)] gives j's preceding vertices at time t - 1, the one-step optimal control, and the weight
        # inspect the reachable set if interested
        if True:
            RS = set()
            for layer in layers[0: self.T + 1]:
                RS = RS | layer
            print('Size of the reachable set is: ', len(RS))
        return predecessors

    def solve(self) -> Tuple[float, List[int], List[int]]:
        """
        Run the solver via dynamic programming. Refer to Algorithm 2 in the paper.


        :return: the optimal value, the optimal control sequence, and the trajectory
        """
        predecessors = self._build_TE_STG()
        M = {}
        U = {}
        F = {}

        def shortest_path(t, j):
            if t == 0:
                return 0
            if (t, j) in F:
                return F[(t, j)]
            i_star = 0
            d_star = float('inf')
            u_i_star = None
            for i, u, w in predecessors[(t, j)]:
                d = shortest_path(t - 1, i) + w
                if d < d_star:
                    i_star = i
                    d_star = d
                    u_i_star = u
            M[(t, j)] = i_star
            F[(t, j)] = d_star
            U[(t-1, i_star, j)] = u_i_star
            return d_star

        T = self.T
        J_star = shortest_path(T + 1, 0)
        u_star = [None] * T
        j = M[(T + 1, 0)]
        t = T
        s = [j]
        while t > 0:
            i = M[(t, j)]
            s.append(i)
            u_star[t - 1] = U[(t - 1, i, j)]
            j = i
            t = t - 1
        s.reverse()
        return J_star, u_star, s






