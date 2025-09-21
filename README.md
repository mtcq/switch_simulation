## Code to accompany: [Can the quantum switch be deterministically simulated?](https://arxiv.org/abs/2409.18202)

#### Marco Túlio Quintino and Jessica Bavaresco

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17170679.svg)](https://doi.org/10.5281/zenodo.17170679)

This is a repository for the article "*Can the quantum switch be deterministically simulated?*, Jessica Bavaresco, Satoshi Yoshida, Tatsuki Odake, Hlér Kristjánsson, Philip Taranto, Mio Murao, and Marco Túlio Quintino, [arXiv:2409.18202 [quant-ph]](https://arxiv.org/abs/2409.18202)".

 The MATLAB code from this repository requires:
- [YALMIP](https://github.com/yalmip/yalmip/) - MATLAB toolbox for optimization modeling
- [MOSEK](https://www.mosek.com) - a software package for solving mathematical optimization problems (under the free personal academic license)
- [SCS](https://www.cvxgrp.org/scs/index.html) - a numerical optimization package for solving large-scale convex cone problems
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory
- [HigherOrderProjectors](https://github.com/mtcq/HigherOrderProjectors) - repository from the paper [arXiv:2305.01247 [quant-ph]](https://arxiv.org/abs/2305.01247) which includes projectors on the linear space spanned by various classes higher-order transformations.

This repository consists of:

#### - ComputerAssisted

- [code_compassist/CertifySDPComputerAssistedProof.m](https://github.com/mtcq/switch_simulation/blob/main/ComputerAssisted/code_compassist/CertifySDPComputerAssistedProof.m): Function that provides a computer-assisted proof Function which constructs a computer assisted proof to certify that the quantum switch cannot be simulated by combs of a given scenario. <br> When dealing with 4-slot combs, this function makes use of Algorithm 3.1. (Veriﬁcation of positive deﬁniteness) from the paper by S.M. Rump, *Verification of Positive Definiteness*, BIT Numerical Mathematics 46, 433–452 (2006), ([https://doi.org/10.1007/s10543-006-0056-1](https://doi.org/10.1007/s10543-006-0056-1)).

- [code_compassist/CertifySDPHighPrecision.m](https://github.com/mtcq/switch_simulation/blob/main/ComputerAssisted/code_compassist/CertifySDPHighPrecision.m): Function that checks with high-precision that the quantum switch cannot be simulated by combs of a given scenario.

- [inputdata_compassist/](https://github.com/mtcq/switch_simulation/tree/main/ComputerAssisted/inputdata_compassist): Folder containing input data required for the computer assisted proof and high-precision code.

- [outputdata_compassist/](https://github.com/mtcq/switch_simulation/tree/main/ComputerAssisted/outputdata_compassist): Folder containing output data obtained with the computer the computer-assisted proof code. That is, these files form a proof certificate for Theorem 2 of [the paper](https://arxiv.org/abs/2409.18202).

#### - SDP

- [linkFAST.m](https://github.com/mtcq/switch_simulation/blob/main/SDP/linkFAST.m): Function that computes the link product of two matrices, $A \in L(\mathbb{C}^{d_{in}})$ and $B \in L(\mathbb{C}^{d_{in}}\otimes\mathbb{C}^{d_{out}})$, returning $C = A*B = tr_{in}((A^T \otimes I) B) \in \mathbb{C}^{d_{out}}$.

- [probabilistic_approximate/](https://github.com/mtcq/switch_simulation/tree/main/SDP/probabilistic_approximate): Folder containing SDP functions for computing the maximum probability of success of an approximate (bounded-error) simulation of the quantum switch in a partly restricted scenario (fixed input systems). Contains functions for the primal problem with simulations using quantum combs and QC-CCs that have 2 and 3 slots. Used to create Fig. 4 of [the paper](https://arxiv.org/abs/2409.18202).

- [probabilistic_exact/](https://github.com/mtcq/switch_simulation/tree/main/SDP/probabilistic_exact): Folder containing SDP functions for computing the maximum probability of success of an exact simulation of the quantum switch in a restricted scenario (fixed input systems, discarded output target system). Contains functions for the primal and dual problems with simulations using quantum combs that have 2, 3 and 4 slots.

#### - data

- [basis/](https://github.com/mtcq/switch_simulation/tree/main/data/basis): Folder containing data files corresponding to bases for the subspace of k identical copies of qubit channels, for $k\in${1,2,3,4}.

- [p=1/](https://github.com/mtcq/switch_simulation/tree/main/data/p%3D1): Folder containing data files corresponding to quantum combs that simulate the action of the quantum switch, in a restricted scenario (fixed input systems, discarded output target system), for identical qubit channels in the order AAAA, and for qubit unitary channels in the order AABB and BAAA.


