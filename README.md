# Finite-Horizon Optimal Control of Boolean Control Networks
This repository contains the code for the paper: Shuhua Gao, Changkai Sun, Cheng Xiang, Kairong Qin, and Tong Heng Lee. *Finite-Horizon Optimal Control of Boolean Control Networks: A Unified Graph-Theoretical Approach* published on *IEEE transactions on neural networks and learning systems*.

**Organization**

The core algorithm implementation is in the folder *src/algorithm*.  We use three examples to show how to use this algorithm in the *src* folder, which correspond to three examples in the paper. The benchmark experiment with the Ara operon network is given in *src/benchmark_ara_operon.py*. The network transition matrix L of this network is presented in *src/ara_operon.txt*.

**Requirement**

Python 3.6 or higher.

**How to run**

+ Download or clone this repository to your local computer.

+ In the command line  (such as *cmd*, *PowerShell* on Windows or *terminal* on Ubuntu), go into the *src* folder

+ To run an example, just type `python ./example_3.py`. (Of course, you can also use certain IDEs like PyCharm.)

**How to cite this work**

Gao, S., Sun, C., Xiang, C., Qin, K., & Lee, T. H. (2020). Finite-Horizon Optimal Control of Boolean Control Networks: A Unified Graph-Theoretical Approach. IEEE transactions on neural networks and learning systems, PP, 10.1109/TNNLS.2020.3027599. Advance online publication. https://doi.org/10.1109/TNNLS.2020.3027599
