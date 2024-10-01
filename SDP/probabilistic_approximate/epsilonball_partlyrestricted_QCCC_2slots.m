% Author: Jessica Bavaresco, jessicabavaresco.github.io
% Requires: QETLAB and all linkFAST.m
% Last update: 01 Oct 2024

function [W,Ws,p] = epsilonball_partlyrestricted_QCCC_2slots(J_A,J_B,d,epsilon)

% Compared to the paper, here we use a relabelling (i,j) -> i of variables Rij -> Ri and {J^A_i, J^B_j} -> {J^A_i, J^B_j}. 

% INPUT: J_A and J_B are N pairs of input channels, with size(J_A) = size(J_B) = [d^2 d^2 N]
%        d the dimension of the target systems
%        epsilon is the fixed approximation parameter   

% OUTPUT: W is the optimal QC-CC simulation
%         Ws is the QC-CC instrument element associated to the success outcome
%         p is the maximum probability of success  


N = size(J_A,3); 
   
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


% Ws_ABF, Ws_BAF, W_ABF, W_BAF in I1 O1 I2 O2 F (here F = to co)
% S in I1 O1 I2 O2 F

Ws_ABF = sdpvar(((d*d)^2)*d*2,((d*d)^2)*d*2,'symmetric','real');
Ws_BAF = sdpvar(((d*d)^2)*d*2,((d*d)^2)*d*2,'symmetric','real');
W_ABF  = sdpvar(((d*d)^2)*d*2,((d*d)^2)*d*2,'symmetric','real');
W_BAF  = sdpvar(((d*d)^2)*d*2,((d*d)^2)*d*2,'symmetric','real');
S      = sdpvar(((d*d)^2)*d*2,((d*d)^2)*d*2,'symmetric','real');

F = [];

for i=1:N
    F = F + [linkFAST(kron(J_A(:,:,i),J_B(:,:,i)),(Ws_ABF+Ws_BAF),[(d*d)^2 d*2])==linkFAST(kron(J_A(:,:,i),J_B(:,:,i)),S,[(d*d)^2 d*2])];
end


F = F + [Ws_ABF>=0,W_ABF-Ws_ABF>=0,W_ABF==ProjSeqSuperChannel(W_ABF,[1 d d d d d*2])];
F = F + [Ws_BAF>=0,W_BAF-Ws_BAF>=0,W_BAF==ProjSeqSuperChannel(W_BAF,[1 d d d d d*2])];
F = F + [trace(W_ABF+W_BAF)==d^2];

F = F + [S>=0,S==ProjGenSuperChannel(S,[1 d d d d d*2])];
F = F + [trace(S*SWITCH)>=(d*d)*trace(S)*(1-epsilon)];

p = trace(S)/(d*d);


solution = solvesdp(F,-p,sdpsettings('solver','mosek','verbose',1,'cachesolvers',1));


p = double(p);
Ws = double(Ws_ABF+Ws_BAF);
W = double(W_ABF+W_BAF);
