"""
Fixed-Destination optimal control

Author: Gao Shuhua
Date: 2019/10/26
"""
from typing import Iterable, Callable, Tuple, List
from operator import itemgetter
from collections import defaultdict
from .bool_net import BooleanNetwork, STG
from .fheap import FibonacciHeap

class STGPlus(STG):
    """
    The extended state transition graph G^+.
    """
    def __init__(self, net: BooleanNetwork, x0: int, Cx: Iterable[int], Cu: Callable[[int], Iterable[int]],
                 Omega: Iterable[int],
                 h: Callable[[int, int], float], g: Callable[[int, int, int], float]
                 ):
        super().__init__(net, x0, Cx, Cu)
        if isinstance(Omega, set):
            self.Omega = Omega
        else:
            self.Omega = set(Omega)
        self.edge_info = {}
        if g is not None:
            self.g = g
        else:
            self.g = lambda x, u, t: 0
        if h is not None:
            self.h = h
        else:
            self.h = lambda x, t: 0
        self._extend()

    def _extend(self):
        """
        (1) introduce the pseudo-state (2) compute the weight and the associated control of each edge
        """
        for i, si in self.adj_list.items():
            for j, us in si.items():
                # an edge i -> j driven by control (a list of multiple possible control inputs)
                ou = min(((k, self.g(i, k, 0)) for k in us), key=itemgetter(1))
                self.edge_info[(i, j)] = ou

        for i in self.Omega:
            if i in self.vertices: # a terminal state
                self.adj_list[i][0] = []    # no control
                self.edge_info[(i, 0)] = (None, self.h(i, 0))
        self.vertices.add(0)


class TimeInvariantCostSolver:
    """
    The solver for fixed-destination optimal control with time-invariant costs, i.e., case 1 of problem 2.
    See Algorithm 3 in the paper.
    """

    def __init__(self, net: BooleanNetwork, x0: int, Cx: Iterable[int], Cu: Callable[[int], Iterable[int]],
                 Omega: Iterable[int],
                 h: Callable[[int, int], float], g: Callable[[int, int, int], float]
                 ):
        self.stgp = STGPlus(net, x0, Cx, Cu, Omega, h, g)   #STG+

    def solve(self) -> Tuple[float, List[int], List[int]]:
        """
        Solve the problem using modified Dijkstra’s algorithm

        :return: the optima value and the optimal control sequence and the state trajectory
        """
        M = {}
        nodes = {}
        Q = FibonacciHeap()
        for i in self.stgp.vertices:
            if i == self.stgp.x0:
                nodes[i] = Q.insert((0, i))
            else:
                nodes[i] = Q.insert((float('inf'), i))
        while Q.total_nodes > 0:
            min_node = Q.extract_min()
            di, i = min_node.data
            if i == 0:
                J_star = di
                break
            for j in self.stgp.adj_list[i]:
                dj = di + self.stgp.edge_info[(i, j)][1]
                if dj < nodes[j].data[0]:
                    Q.decrease_key(nodes[j], (dj, j))
                    M[j] = i
        # reconstruct the optimal control sequence
        u_star = []
        j = 0
        s = []
        while j != self.stgp.x0:
            i = M[j]
            s.append(i)
            if j != 0:
                u_star.append(self.stgp.edge_info[(i, j)][0])
            j = i
        u_star.reverse()
        s.reverse()
        return J_star, u_star, s


class TEDSTG:
    """
    The TED-STG.
    """
    def __init__(self, net: BooleanNetwork, x0: int, Cx: Iterable[int], Cu: Callable[[int], Iterable[int]],
                 Omega: Iterable[int],
                 h: Callable[[int, int], float], g: Callable[[int, int, int], float]
                 ):
        self.net = net
        self.x0 = x0
        self.Cx = Cx
        self.Cu = Cu
        if isinstance(Omega, set):
            self.Omega = Omega
        else:
            self.Omega = set(Omega)
        self.edge_info = {}
        if g is not None:
            self.g = g
        else:
            self.g = lambda x, u, t: 0
        if h is not None:
            self.h = h
        else:
            self.h = lambda x, t: 0
        self.edge_info = {}
        self.adj_list = defaultdict(list)
        self.layers = []
        self._build()

    def _build(self):
        stg = STG(self.net, self.x0, self.Cx, self.Cu)
        Z = len(stg.vertices)
        layers = [set() for _ in range(Z + 1)]
        layers[0].add(self.x0)
        for t in range(Z - 1):
            for i in layers[t]:
                for j, us in stg.adj_list[i].items():
                    self.adj_list[(t, i)].append(j)
                    # optimal control from i to j
                    u, w = min(((k, self.g(i, k, t)) for k in us), key=itemgetter(1))
                    self.edge_info[(t, i, j)] = (u, w)
                    layers[t + 1].add(j)
        layers[Z].add(0)
        # handle V_Z and E0
        for t in range(Z):
            for i in self.Omega:
                if i in layers[t]:
                    # connect i to \delta^0
                    self.adj_list[(t, i)].append(0)
                    # no control, and the weight is given by the terminal cost
                    self.edge_info[(t, i, 0)] = (None, self.h(i, t))
        self.layers = layers



class TimeVariantCostSolver:
    """
    The solver for fixed-destination optimal control with time-variant costs, i.e., case 1 of problem 2.
    See Algorithm 3 in the paper.
    """

    def __init__(self, net: BooleanNetwork, x0: int, Cx: Iterable[int], Cu: Callable[[int], Iterable[int]],
                 Omega: Iterable[int],
                 h: Callable[[int, int], float], g: Callable[[int, int, int], float]
                 ):
        self.tedstg = TEDSTG(net, x0, Cx, Cu, Omega, h, g)

    def solve(self) -> Tuple[float, List[int], List[int]]:
        """
        Solve the problem using modified Dijkstra’s algorithm

        :return: the optima value and the optimal control sequence and the state trajectory
        """
        M = {}
        nodes = {}
        Q = FibonacciHeap()
        Z = len(self.tedstg.layers) - 1
        for t in range(len(self.tedstg.layers)):
            for i in self.tedstg.layers[t]:
                if t == 0:
                    nodes[(t, i)] = Q.insert((0, t, i))
                else:
                    nodes[(t, i)] = Q.insert((float('inf'), t, i))
        while Q.total_nodes > 0:
            min_node = Q.extract_min()
            di, t, i = min_node.data
            if i == 0:
                J_star = di
                break
            for j in self.tedstg.adj_list[(t, i)]:
                d = di + self.tedstg.edge_info[(t, i, j)][1]
                tj = Z if j == 0 else t + 1  # the pseudo-state is represented by \delta_{N, Z}^0 here
                node = nodes[(tj, j)]
                if d < node.data[0]:
                    Q.decrease_key(node, (d, tj, j))
                    M[(tj, j)] = (t, i)
        # reconstruct the control sequence
        u_star = []
        tj = Z
        j = 0
        s = []
        while j != self.tedstg.x0:
            ti, i = M[(tj, j)]
            s.append(i)
            if j != 0:  # the pseudo-state needs no control
                u_star.append(self.tedstg.edge_info[(ti, i, j)][0])
            j = i
            tj = ti
        u_star.reverse()
        s.reverse()
        return J_star, u_star, s
