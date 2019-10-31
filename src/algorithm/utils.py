"""
Some utility functions

Author: Gao Shuhua
Date: 2019/10/25
"""
from typing import Tuple, List


def read_network(file: str) -> Tuple[int, int, List[int]]:
    """
    Read a Boolean network from a text file:
        Line 1: number of state variables
        Line 2: number of control inputs
        Line 4: transition matrix of the network (linear representation of a logical matrix)
        ...

    :param file: a text file
    :return: (n, m, L), where
        n: number of state variables
        m: number of control inputs
        Ls: a list of transition matrices, each for one sub-system
    """
    with open(file, 'r') as f:
        n = int(f.readline().strip())
        m = int(f.readline().strip())
        N = 2 ** n
        M = 2 ** m
        line = f.readline().strip()
        assert line, f'{w} transition matrices must be provided!'
        numbers = line.split()
        assert len(numbers) == M * N, f'The transition matrix must have {M * N} columns'
        L = [int(num) for num in numbers]
        return n, m, L

