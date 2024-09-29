## Code to accompany: [Can the quantum switch be deterministically simulated?](https://arxiv.org/abs/xxxx.xxxx)

#### Marco Túlio Quintino and Jessica Bavaresco


This is a repository for the article [Can the quantum switch be deterministically simulated?](https://arxiv.org/abs/xxxx.xxxx), from Jessica Bavaresco, Satoshi Yoshida, Tatsuki Odake, Hlér Kristjánsson, Philip Taranto, Mio Murao, and Marco Túlio Quintino.

 The MATLAB code from this repository requires:
- [YALMIP](https://github.com/yalmip/yalmip/) - MATLAB toolbox for optimization modeling
- [QETLAB](http://www.qetlab.com/) - a free MATLAB toolbox for quantum entanglement theory
- [HigherOrderProjectors](https://github.com/mtcq/HigherOrderProjectors) - repository from the article [https://arxiv.org/abs/2305.01247](https://arxiv.org/abs/2305.01247) which includdes projectors on the linear space spanned by various classes higher-order transformations.

This repository consists of:

- [CertifySDPComputerAssistedProof.m](https://github.com/mtcq/switch_simulation/blob/main/code_compassist/CertifySDPComputerAssistedProof.m): Function that provides a computer-assisted proof Function which constructs a computer assisted proof to certify that the quantum switch cannot be simulated by combs of a given scenario. <br> When dealing with 4-slot combs, this function makes use of Algorithm 3.1. (Veriﬁcation of positive deﬁniteness) from the article Verification of Positive Definiteness, written by S.M. Rump ([https://doi.org/10.1007/s10543-006-0056-1](https://doi.org/10.1007/s10543-006-0056-1)).

- [CertifySDPHighPrecision.m](https://github.com/mtcq/switch_simulation/blob/main/code_compassist/CertifySDPHighPrecision.m): Function which checks with high-precision that the quantum switch cannot be simulated by combs of a given scenario.

-  [inputdata_compassist](https://github.com/mtcq/switch_simulation/tree/main/inputdata_compassist): Folder which contains input data required for the computer assisted proof and high-precision code.

- [outputdata_compassist.m](https://github.com/mtcq/switch_simulation/tree/main/outputdata_compassist): Folder which contains output data obtained with the computer the computer-assisted proof code. That is, these files form a proof certificate for Theorem 2 of the article [Can the quantum switch be deterministically simulated?](https://arxiv.org/abs/xxxx.xxxx).

