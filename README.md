## Code to accompany: [Can the quantum switch be deterministically simulated?](https://arxiv.org/abs/2409.18202)

#### Marco Túlio Quintino and Jessica Bavaresco

This is a repository for the article "*Can the quantum switch be deterministically simulated?*, Jessica Bavaresco, Satoshi Yoshida, Tatsuki Odake, Hlér Kristjánsson, Philip Taranto, Mio Murao, and Marco Túlio Quintino, [arXiv:2409.18202 [quant-ph]](https://arxiv.org/abs/2409.18202)".

 The MATLAB code from this repository requires:
- [YALMIP](https://github.com/yalmip/yalmip/) - MATLAB toolbox for optimization modeling
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory
- [HigherOrderProjectors](https://github.com/mtcq/HigherOrderProjectors) - repository from the article  [arXiv:2305.01247 [quant-ph]](https://arxiv.org/abs/2305.01247) which includes projectors on the linear space spanned by various classes higher-order transformations.

This repository consists of:

#### - ComputerAssisted

- [code_compassist/CertifySDPComputerAssistedProof.m](https://github.com/mtcq/switch_simulation/blob/main/ComputerAssisted/code_compassist/CertifySDPComputerAssistedProof.m): Function that provides a computer-assisted proof Function which constructs a computer assisted proof to certify that the quantum switch cannot be simulated by combs of a given scenario. <br> When dealing with 4-slot combs, this function makes use of Algorithm 3.1. (Veriﬁcation of positive deﬁniteness) from the article Verification of Positive Definiteness, written by S.M. Rump ([https://doi.org/10.1007/s10543-006-0056-1](https://doi.org/10.1007/s10543-006-0056-1)).

- [code_compassist/CertifySDPHighPrecision.m](https://github.com/mtcq/switch_simulation/blob/main/ComputerAssisted/code_compassist/CertifySDPHighPrecision.m): Function which checks with high-precision that the quantum switch cannot be simulated by combs of a given scenario.

- [inputdata_compassist/](https://github.com/mtcq/switch_simulation/tree/main/ComputerAssisted/inputdata_compassist): Folder which contains input data required for the computer assisted proof and high-precision code.

- [outputdata_compassist/](https://github.com/mtcq/switch_simulation/tree/main/ComputerAssisted/outputdata_compassist): Folder which contains output data obtained with the computer the computer-assisted proof code. That is, these files form a proof certificate for Theorem 2 of the article [Can the quantum switch be deterministically simulated?](https://arxiv.org/abs/2409.18202).

#### - SDP

- [linkFAST.m](https://github.com/mtcq/switch_simulation/tree/main/SDP/linkFAST.m): Function that computes the link product of two matrices, $A \in I$ and $B \in I\otimes O$, returning $C = A*B \in O$.

- 

### - data

- [basis/](https://github.com/mtcq/switch_simulation/tree/main/data/basis): Folder which contains data files corresponding to bases for the subspace of k identical copies of qubit channels, for $k\in\{1,2,3,4\}.$

- [p=1/](https://github.com/mtcq/switch_simulation/tree/main/data/p%3D1): Folder which contains data files corresponding to quantum combs that simulate the action of the quantum switch, in a restricted scenario (fixed input systems, discarded output target system), for identical qubit channels in the order AAAA, and for qubit unitary channels in the order AABB and BAAA.


