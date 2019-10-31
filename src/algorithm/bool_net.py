"""
A Boolean network

Author: Gao Shuhua
Date: 2019/10/24
"""
from typing import Callable, List, Set, Iterable, Dict
from collections import defaultdict
from collections import deque


class BooleanNetwork:
    """
    A Boolean network
    """
    def __init__(self, L: List[int], m: int, n: int):
        """
        Initialize a Boolean network.

        :param L: network transition matrix in a linear form
        :param m: number of control inputs
        :param n: number of state variables
        """
        M = 2 ** m
        N = 2 ** n
        assert len(L) == M * N
        self.L = L
        self.m = m
        self.n = n
        self.M = M
        self.N = N

    def successors(self, i: int, Cx: Iterable[int]=None, Cu: Callable[[int], Iterable[int]]=None) -> Dict[int, List[int]]:
        """
        Get the succeeding states of `i` under constraints.

        :param i: current state
        :param Cx: state constraints, only states in Cx are allowed
        :param Cu: control constraints, Cu(i) gives the control that are allowed at i
        :return: successors and the control inputs that attain each successor
        """
        successors = defaultdict(list)
        if Cu is None:  # no input constraints at all
            us = (k for k in range(1, self.M + 1, 1))
        else:
            us = Cu(i)
        for k in us:
            blk = self.L[(k - 1) * self.N: k * self.N]
            j = blk[i - 1]
            if Cx is None or j in Cx:
                successors[j].append(k)
        assert len(successors) > 0, 'Invalid constraint settings! Each admissible state must have successors.'
        return successors


class STG:
    """
    State transition graph.
    Only the topology of the directed graph is built, and no weights are assigned.
    See Algorithm 1 in the paper.
    """
    def __init__(self, net: BooleanNetwork, x0: int, Cx: Iterable[int], Cu: Callable[[int], Iterable[int]]):
        self.net = net
        self.x0 = x0
        self.Cx = Cx
        self.Cu = Cu
        self.adj_list = {}  # adj_list[i] is a dict, where j -> us means j is a successor of i attained with controls us
        self.vertices = None
        self._build()

    def _build(self):
        Q = deque([self.x0])
        R = set()
        while Q:
            i = Q.popleft()
            self.adj_list[i] = self.net.successors(i, self.Cx, self.Cu)
            for j in self.adj_list[i]:
                if j not in R:
                    R.add(j)
                    Q.append(j)
        self.vertices = R



