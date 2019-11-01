"""

Author: Gao Shuhua
Date: 2019/10/28
"""
from typing import List, Callable
import copy
from operator import itemgetter
from .fheap import FibonacciHeap
import numpy as np
import time


class AbstractSolver:
    """
    A base class to provide a common interface fro various methods.
    This is used for the benchmark purpose: only time-invariant stage cost and zero terminal cost are supported. Besides,
    there are no constraints.
    """
    _INF = float('inf')

    def __init__(self, L: List[int], m: int, n: int, g: Callable[[int, int], float]):
        self.L = L
        self.m = m
        self.n = n
        self.M = 2 ** m
        self.N = 2 ** n
        self.g = g

    def solve(self, x0: int, xd: int, T: int=None) -> float:
        """
        Solve and return the optimal value

        :param x0: initial state
        :param xd: terminal state
        :param T: if None, then a fixed-destination problem; otherwise, a fixed-time problem.
        :return: optimal value
        """
        raise NotImplementedError()


class LiFangfeiSolver(AbstractSolver):
    """
    [1] F. Li and X. Lu, “Minimum energy control and optimal-satisfactory
    control of boolean control network,” Physics Letters A, vol. 377, no. 43,
    pp. 3112–3118, 2013.
    """
    def __init__(self, L: List[int], m: int, n: int, g: Callable[[int, int], float]):
        super().__init__(L, m, n, g)
        self.D = self._compute_cost_matrix()    # one-step cost matrix

    def _compute_cost_matrix(self):
        D = [[float('inf')] * (self.N + 1) for _ in range(self.N + 1)]
        for i in range(1, self.N + 1):
            for k in range(1, self.M + 1):
                # i --> j with control k
                blk = self.L[(k - 1) * self.N: k * self.N]
                j = blk[i - 1]
                D[i][j] = min(D[i][j], self.g(i, k))
        return D

    def solve(self, x0: int, xd: int, T: int=None) -> float:
        """
        Solve min energy control with a given horizon length or an unknown length if `T` is None.

        :param x0: initial state
        :param xd: destination state
        :param T: horizon length or None
        :return: optimal value
        """
        D = copy.deepcopy(self.D)
        next_D = copy.deepcopy(D)

        if T is None:
            # unspecified horizon length: Algorithm 3.1 in [1]
            for d in range(self.N - 1):
                for i in range(1, self.N + 1):
                    for j in range(1, self.N + 1):
                        next_D[i][j] = min(D[i][j], min(D[i][alpha] + self.D[alpha][j] for alpha in range(1, self.N + 1)))
                D, next_D = next_D, D
        else:
            # fixed horizon length: Algorithm 3.2 in [1]
            for d in range(T - 1):
                for i in range(1, self.N + 1):
                    for j in range(1, self.N + 1):
                        next_D[i][j] = min(D[i][alpha] + self.D[alpha][j] for alpha in range(1, self.N + 1))
                D, next_D = next_D, D

        return D[x0][xd]


class EttoreFornasiniSolver(AbstractSolver):
    """
     E. Fornasini and M. E. Valcher, “Optimal control of boolean control
    networks,” IEEE Transactions on Automatic Control, vol. 59, no. 5, pp.
    1258–1270, 2013.
    """

    def __init__(self, L: List[int], m: int, n: int, g: Callable[[int, int], float]):
        super().__init__(L, m, n, g)

    def solve(self, x0: int, xd: int, T: int=None) -> float:
        assert isinstance(T, int), 'ONLY horizon of a fixed length is supported'
        M = self.M
        N = self.N
        # compute the stage cost vector c (written in a matrix form for convenience)
        c = [[None] * (N + 1) for _ in range(M + 1)]
        for k in range(1, M + 1):
            for i in range(1, N + 1):
                c[k][i] = self.g(i, k)
        # only xd has zero terminal cost; others all infinity
        cf = [self._INF] * (N + 1)
        cf[xd] = 0
        # main loop
        m = self._compute_m_without_optimization(cf, c, T, M, N, self.L)
        return m[0][x0]

    def _compute_m_without_optimization(self, cf, c, T, M, N, L):
        """
        Calculate the vector m in the recursive algorithm.
        Note the cost vectors cf and c start indexing from 1.
        Refer to the algorithm given on page 1261.

        :return: a 2d list
        """
        assert len(cf) == N + 1
        assert len(c) == M + 1 and len(c[1]) == N + 1
        assert len(L) == M * N
        # initialization
        m = [[None] * (N + 1) for _ in range(T)]
        m.append(cf)
        # recursion
        for t in range(T - 1, -1, -1):
            for j in range(1, N + 1):
                values = []
                for i in range(1, M + 1):
                    Li = L[(i - 1) * N: i * N]  # the ith-block
                    values.append(c[i][j] + m[t + 1][Li[j - 1]])
                m[t][j] = min(values)
        return m


class ZhuQunxiSolver(AbstractSolver):
    """
     Q. Zhu, Y. Liu, J. Lu, and J. Cao, “On the optimal control of boolean
    control networks,” SIAM Journal on Control and Optimization, vol. 56,
    no. 2, pp. 1321–1341, 2018.
    """

    def __init__(self, L: List[int], m: int, n: int, g: Callable[[int, int], float]):
        super().__init__(L, m, n, g)
        self.D = [[self._INF] * (self.N + 1) for _ in range(self.N + 1)]
        self.C = [[None] * (self.N + 1) for _ in range(self.N + 1)]
        for i in range(1, self.N + 1):
            for k in range(1, self.M + 1):
                # i --> j with control k
                blk = self.L[(k - 1) * self.N: k * self.N]
                j = blk[i - 1]
                cik = self.g(i, k)
                if cik < self.D[i][j]:
                    self.D[i][j] = cik
                    self.C[i][j] = k

    def _binary_decompose(self, T: int) -> List[int]:
        assert T > 0
        s = []
        while T > 0:
            s.append(T & 1)
            T = T >> 1
        return s

    def _circle(self, A, B):
        """
        The circle operator.

        :param A: n-by-n matrix
        :param B: n-by-n matrix
        :return: a n-by-n matrix
        """
        assert len(A) == len(B)
        n = len(A) - 1
        C = copy.deepcopy(A)
        for i in range(1, n + 1):
            for j in range(1, n + 1):
                C[i][j] = min(A[i][k] + B[k][j] for k in range(1, n + 1))
        return C

    def _stp(self, i, p, j, q):
        """
        STP of \delta_p^i and \delta_q^j

        :return: t of \delta_{p*q}^t
        """
        if i is None or j is None:
            return None
        return (i - 1) * q + j


    def solve(self, x0: int, xd: int, T: int=None) -> float:
        assert isinstance(T, int), 'ONLY horizon of a fixed length is supported'
        M = self.M
        N = self.N
        betas = self._binary_decompose(T)
        l = len(betas) - 1
        ts = 0
        C = {}
        C[1] = self.C
        D = {}
        D[1] = self.D
        old_ts = 0   # t_{s-1}
        ts = 0   # t_s
        for s in range(l + 1):
            # step 1
            ts = old_ts + 2 ** s * betas[s]
            # step 2
            p = 2 ** s
            q = 2 ** (s + 1)
            if betas[s] == 1 and ts not in D:
                D[ts] = self._circle(D[old_ts], D[p])
                C[ts] = [[None] * (self.N + 1) for _ in range(self.N + 1)]
                for i in range(1, N + 1):
                    for j in range(1, N + 1):
                        k_star = min(((k, D[old_ts][i][k] + D[p][k][j]) for k in range(1, N + 1)), key=itemgetter(1))[0]
                        if D[ts][i][j] == self._INF:
                            C[ts][i][j] = None
                        else:
                            C[ts][i][j] = self._stp(C[old_ts][i][k_star], M ** old_ts, C[p][k_star][j], M ** p)
            # step 3
            if q <= T:
                D[q] = self._circle(D[p], D[p])
                C[q] = [[None] * (self.N + 1) for _ in range(self.N + 1)]
                for i in range(1, N + 1):
                    for j in range(1, N + 1):
                        k_star = min(((k, D[p][i][k] + D[p][k][j]) for k in range(1, N + 1)), key=itemgetter(1))[0]
                        if D[q][i][j] == self._INF:
                            C[q][i][j] = None
                        else:
                            C[q][i][j] = self._stp(C[p][i][k_star], M ** p, C[p][k_star][j], M ** p)
            old_ts = ts
        # we only allow a single destination state
        return D[T][x0][xd]


class CuiXingbangSolver(AbstractSolver):
    """
     X. Cui, J.-E. Feng, and S. Wang, “Optimal control problem of boolean
    control networks: A graph-theoretical approach,” in 2018 Chinese Con-
    trol And Decision Conference (CCDC). IEEE, 2018, pp. 4511–4516.

    The above paper adopts a graph similar to our TET-STG, but it doesn't consider time-variant costs. Besides, it uses
    Dijkstra's algorithm to find the shortest path instead of our more efficient dynamic programming method.

    This implementation here only solves the minimum-energy problem. As for the minimum-time one, the method in that
    paper is essentially the same as ours.
    """
    def __init__(self, L: List[int], m: int, n: int, g: Callable[[int, int], float]):
        super().__init__(L, m, n, g)
        self._successors = {}

    def _get_successors(self, i):
        s = {}
        for k in range(1, self.M + 1):
            blk = self.L[(k - 1) * self.N: k * self.N]
            j = blk[i - 1]
            c = self.g(i, k)
            if c < s.get(j, self._INF):
                s[j] = c
        return s

    def _build_graph(self, x0, xd, T):
        layers = [set() for _ in range(T + 1)]
        layers[0].add(x0)
        for t in range(1, T + 1):
            for i in layers[t - 1]:
                # s = self._successors(i)
                s = self._get_successors(i)
                self._successors[(t - 1, i)] = s
                layers[t] = layers[t] | s.keys()
        # the final layers only contains the single destination state
        assert xd in layers[T], 'The given destination state cannot be reached.'
        layers[T] = {xd}
        return layers

    def solve(self, x0: int, xd: int, T: int=None) -> float:
        # Dijkstra's algorithm
        assert isinstance(T, int), 'ONLY horizon of a fixed length is supported'
        layers = self._build_graph(x0, xd, T)
        nodes = {}
        Q = FibonacciHeap()
        for t in range(T + 1):
            for i in layers[t]:
                if t == 0:
                    nodes[(t, i)] = Q.insert((0, t, i))
                else:
                    nodes[(t, i)] = Q.insert((self._INF, t, i))
        while Q.total_nodes > 0:
            min_node = Q.extract_min()
            di, t, i = min_node.data
            if t == T:
                J_star = di
                break
            # update the value of its sucessors
            s = self._successors[(t, i)]
            if t == T - 1:
                if xd in s:
                    d = di + s[xd]
                    node = nodes[(T, xd)]
                    if d < node.data[0]:
                        Q.decrease_key(node, (d, T, xd))
            else:
                for j in s:
                    d = di + s[j]
                    node = nodes[(t + 1, j)]
                    if d < node.data[0]:
                        Q.decrease_key(node, (d, t + 1, j))
        return J_star


class DmitriyLaschovSolver(AbstractSolver):
    """
     D. Laschov and M. Margaliot, “Minimum-time control of boolean
    networks,” SIAM Journal on Control and Optimization, vol. 51, no. 4,
    pp. 2869–2892, 2013.
    """
    def __init__(self, L: List[int], m: int, n: int, g: Callable[[int, int], float]):
        super().__init__(L, m, n, g)

    def _boolean_product(self, P: np.ndarray, Q: np.ndarray) -> np.ndarray:
        m, n = P.shape
        p, q = Q.shape
        assert n == p
        R = np.empty((m, q), np.int8)
        for i in range(m):
            for j in range(q):
                r = False
                for k in range(n):
                    r = r or (P[i][k] and Q[k][j])
                    if r:
                        break
                R[i][j]= r
        return R

    def _boolean_stp(self, P: np.ndarray, Q: np.ndarray):
        m, n = P.shape
        p, q = Q.shape
        alpha = np.lcm(n, p)
        return self._boolean_product(np.kron(P, np.eye(alpha // n)), np.kron(Q, np.eye(alpha // p)))

    def solve(self, x0: int, xd: int, T: int=None):
        # Refer to Theorem 4 therein.
        assert T is None, 'ONLY minimum-time problem is supported'
        M = self.M
        N = self.N
        # change L to its matrix form (N-by-MN)
        Lm = np.zeros((N, M * N), dtype=np.int8)
        for j, i in enumerate(self.L):
            Lm[i - 1][j] = 1
        Q = self._boolean_stp(Lm, np.ones((M, 1), dtype=np.int8))
        # change x0 and xd to its vector form
        z = np.zeros((N, 1), np.int8)
        z[xd - 1] = 1
        temp = x0
        x0 = np.zeros((N, 1), np.int8)
        x0[temp - 1] = 1
        # main loop
        eta = z
        N_star = -1
        for k in range(N):
            r = (eta.T @ x0).item()
            if r == 1:
                N_star = k
                break
            eta = self._boolean_product(Q.T, eta)
        return N_star









