% Author: Jessica Bavaresco, jessicabavaresco.github.io
% Requires: QETLAB and all linkFAST.m
% Last update: 01 Oct 2024

function [W,Ws,p] =  epsilonball_partlyrestricted_QCCC_3slots(J_A,J_B,d,epsilon)

% Compared to the paper, here we use a relabelling (i,j) -> i of variables Rij -> Ri and {J^A_i, J^B_j} -> {J^A_i, J^B_j}. 

% INPUT: J_A and J_B are N pairs of input channels, with size(J_A) = size(J_B) = [d^2 d^2 N]
%        d the dimension of the target systems
%        epsilon is the fixed approximation parameter   

% OUTPUT: W is the optimal QC-CC simulation
%         Ws is the QC-CC instrument element associated to the success outcome
%         p is the maximum probability of success  


N = size(J_A,3); 

dAi = d;
dAo = d;
dBi = d;
dBo = d;
dCi = d;
dCo = d;
dF  = 2*d;
        
alpha = 1/sqrt(2);
beta  = sqrt(1-alpha^2);
% control input system in state |+> in dimension 2

target = zeros(d,1);
target(1,1) = 1;
% target input system in state |0> in dimension d

phi = sqrt(d)*MaxEntangled(d);
% Choi vector of identity channel in dimension d

SWITCH = alpha*Tensor(target,phi,phi,[1;0]) + beta*PermuteSystems(Tensor(target,phi,phi,[0;1]),[3 4 1 2 5 6],[d d d d d 2]);
% SWITCH in Ai Ao Bi Bo to co 
% dim(Ai)=dim(Ao)=dim(Bi)=dim(Bo)=dim(to)=d, dim(co)=2

SWITCH = SWITCH*SWITCH';
% partly restricted simulation scenario


yalmip('clear')


% all W and Ws in I1 O1 I2 O2 I3 O3 F (here F = to co)
% S in I1 O1 I2 O2 F

W_ABCF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
W_ACBF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
W_BACF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
W_BCAF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
W_CABF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
W_CBAF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');

Ws_ABCF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
Ws_ACBF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
Ws_BACF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
Ws_BCAF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
Ws_CABF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');
Ws_CBAF = sdpvar(((d*d)^3)*d*2,((d*d)^3)*d*2,'symmetric','real');

S = sdpvar(((d*d)^2)*d*2,((d*d)^2)*d*2,'symmetric','real');

F = [];

for i=1:N
    F = F + [linkFAST(Tensor(J_A(:,:,i),J_B(:,:,i),J_A(:,:,i)),(Ws_ABCF+Ws_ACBF+Ws_BACF+Ws_BCAF+Ws_CABF+Ws_CBAF),[(d*d)^3 d*2])==linkFAST(kron(J_A(:,:,i),J_B(:,:,i)),S,[(d*d)^2 d*2])];
end

F = F + [trace(W_ABCF+W_ACBF+W_BACF+W_BCAF+W_CABF+W_CBAF)==d^3];

F = F + [Ws_ABCF>=0,W_ABCF-Ws_ABCF>=0];
F = F + [Ws_ACBF>=0,W_ACBF-Ws_ACBF>=0];
F = F + [Ws_BACF>=0,W_BACF-Ws_BACF>=0];
F = F + [Ws_BCAF>=0,W_BCAF-Ws_BCAF>=0];
F = F + [Ws_CABF>=0,W_CABF-Ws_CABF>=0];
F = F + [Ws_CBAF>=0,W_CBAF-Ws_CBAF>=0];

% constraints for W_AF = W_ABCF + W_ACBF
F = F + [PartialTrace(W_ABCF+W_ACBF,[2],[(dAi*dAo) (dBi*dBo)*(dCi*dCo)*dF])==kron(PartialTrace(W_ABCF+W_ACBF,[2],[dAi dAo*(dBi*dBo)*(dCi*dCo)*dF]),eye(dAo)/dAo)];

% constraints for W_BF = W_BACF + W_BCAF
F = F + [PartialTrace(W_BACF+W_BCAF,[1 3],[(dAi*dAo) (dBi*dBo) (dCi*dCo)*dF])==kron(PartialTrace(W_BACF+W_BCAF,[1 3],[(dAi*dAo) dBi dBo*(dCi*dCo)*dF]),eye(dBo)/dBo)];

% constraints for W_CF = W_CABF + W_CBAF
F = F + [PartialTrace(W_CABF+W_CBAF,[1 3],[(dAi*dAo)*(dBi*dBo) (dCi*dCo) dF])==kron(PartialTrace(W_CABF+W_CBAF,[1 3],[(dAi*dAo)*(dBi*dBo) dCi dCo*dF]),eye(dCo)/dCo)];

% constraints for W_ABCF
F = F + [PartialTrace(W_ABCF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF])==kron(PartialTrace(W_ABCF,[2],[(dAi*dAo)*(dBi*dBo)*dCi dCo*dF]),eye(dCo)/dCo)];
F = F + [PartialTrace(W_ABCF,[2],[(dAi*dAo)*(dBi*dBo) (dCi*dCo)*dF])==kron(PartialTrace(W_ABCF,[2],[(dAi*dAo)*dBi dBo*(dCi*dCo)*dF]),eye(dBo)/dBo)];
% for no dynamical order
%F = F + [PartialTrace(W_ABCF,[2],[(dAi*dAo) (dBi*dBo)*(dCi*dCo)*dF])==kron(PartialTrace(W_ABCF,[2],[dAi dAo*(dBi*dBo)*(dCi*dCo)*dF]),eye(dAo)/dAo)];

% constraints for W_ACBF
F = F + [PartialTrace(W_ACBF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF])==TR(PartialTrace(W_ACBF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF]),[2],[(dAi*dAo)*dBi dBo (dCi*dCo)])]; 
F = F + [PartialTrace(W_ACBF,[2 4],[(dAi*dAo) (dBi*dBo) (dCi*dCo) dF])==kron(PartialTrace(W_ACBF,[2 4],[(dAi*dAo) (dBi*dBo) dCi dCo*dF]),eye(dCo)/dCo)];
% for no dynamical order
%F = F + [PartialTrace(W_ACBF,[2],[(dAi*dAo) (dBi*dBo)*(dCi*dCo)*dF])==kron(PartialTrace(W_ACBF,[2],[dAi dAo*(dBi*dBo)*(dCi*dCo)*dF]),eye(dAo)/dAo)];

% constraints for W_BACF
F = F + [PartialTrace(W_BACF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF])==kron(PartialTrace(W_BACF,[2],[(dAi*dAo)*(dBi*dBo)*dCi dCo*dF]),eye(dCo)/dCo)];
F = F + [PartialTrace(W_BACF,[2],[(dAi*dAo)*(dBi*dBo) (dCi*dCo)*dF])==TR(PartialTrace(W_BACF,[2],[(dAi*dAo)*(dBi*dBo) (dCi*dCo)*dF]),[2],[dAi dAo (dBi*dBo)])];
% for no dynamical order
%F = F + [PartialTrace(W_BACF,[1 3],[(dAi*dAo) (dBi*dBo) (dCi*dCo)*dF])==kron(PartialTrace(W_BACF,[1 3],[(dAi*dAo) dBi dBo*(dCi*dCo)*dF]),eye(dBo)/dBo)];

% constraints for W_BCAF
F = F + [PartialTrace(W_BCAF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF])==TR(PartialTrace(W_BCAF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF]),[2],[dAi dAo (dBi*dBo)*(dCi*dCo)])]; 
F = F + [PartialTrace(W_BCAF,[1 3],[(dAi*dAo) (dBi*dBo)*(dCi*dCo) dF])==kron(PartialTrace(W_BCAF,[1 3],[(dAi*dAo) (dBi*dBo)*dCi dCo*dF]),eye(dCo)/dCo)];
% for no dynamical order
%F = F + [PartialTrace(W_BCAF,[1 3],[(dAi*dAo) (dBi*dBo) (dCi*dCo)*dF])==kron(PartialTrace(W_BCAF,[1 3],[(dAi*dAo) dBi dBo*(dCi*dCo)*dF]),eye(dBo)/dBo)];

% constraints for W_CABF
F = F + [PartialTrace(W_CABF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF])==TR(PartialTrace(W_CABF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF]),[2],[(dAi*dAo)*dBi dBo (dCi*dCo)])]; 
F = F + [PartialTrace(W_CABF,[2 4],[(dAi*dAo) (dBi*dBo) (dCi*dCo) dF])==TR(PartialTrace(W_CABF,[2 4],[(dAi*dAo) (dBi*dBo) (dCi*dCo) dF]),[2],[dAi dAo (dCi*dCo)])];
% for no dynamical order
%F = F + [PartialTrace(W_CABF,[1 3],[(dAi*dAo)*(dBi*dBo) (dCi*dCo) dF])==kron(PartialTrace(W_CABF,[1 3],[(dAi*dAo)*(dBi*dBo) dCi dCo*dF]),eye(dCo)/dCo)];

% constraints for W_CBAF
F = F + [PartialTrace(W_CBAF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF])==TR(PartialTrace(W_CBAF,[2],[(dAi*dAo)*(dBi*dBo)*(dCi*dCo) dF]),[2],[dAi dAo (dBi*dBo)*(dCi*dCo)])]; 
F = F + [PartialTrace(W_CBAF,[1 3],[(dAi*dAo) (dBi*dBo)*(dCi*dCo) dF])==TR(PartialTrace(W_CBAF,[1 3],[(dAi*dAo) (dBi*dBo)*(dCi*dCo) dF]),[2],[dBi dBo (dCi*dCo)])]; 
% for no dynamical order
%F = F + [PartialTrace(W_CBAF,[1 3],[(dAi*dAo)*(dBi*dBo) (dCi*dCo) dF])==kron(PartialTrace(W_CBAF,[1 3],[(dAi*dAo)*(dBi*dBo) dCi dCo*dF]),eye(dCo)/dCo)]; 


F = F + [S>=0,S==ProjGenSuperChannel(S,[1 d d d d d*2])];
F = F + [trace(S*SWITCH)>=(d*d)*trace(S)*(1-epsilon)];

p = trace(S)/(d*d);


solution = solvesdp(F,-p,sdpsettings('solver','scs','verbose',1,'cachesolvers',1));
%solution = solvesdp(F,-p,sdpsettings('solver','scs','verbose',1,'cachesolvers',1,'scs.max_iters',10000,'scs.eps_abs',1e-6,'scs.eps_rel',1e-6));


p = double(p);
Ws = double(Ws_ABCF+Ws_ACBF+Ws_BACF+Ws_BCAF+Ws_CABF+Ws_CBAF);
W = double(W_ABCF+W_ACBF+W_BACF+W_BCAF+W_CABF+W_CBAF);
