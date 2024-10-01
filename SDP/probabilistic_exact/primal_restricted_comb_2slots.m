% Author: Jessica Bavaresco, jessicabavaresco.github.io
% Requires: QETLAB and all linkFAST.m
% Last update: 01 Oct 2024

function [C,Cs,p] = primal_restricted_comb_2slots(J_A,J_B,d,order)

% Compared to the paper, here we use a relabelling (i,j) -> i of variables Rij -> Ri and {J^A_i, J^B_j} -> {J^A_i, J^B_j}. 

% INPUT: J_A and J_B are N pairs of input channels, with size(J_A) = size(J_B) = [d^2 d^2 N]
%        d the dimension of the target systems
%        order is a vector with the order in which the comb will act on the input channels, e.g., for the order AB, order = [1 2]

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

SWITCH = PartialTrace(SWITCH*SWITCH',2,[d^4 d 2]);
% SWITCH in Ai Ao Bi Bo co 
% dim(Ai)=dim(Ao)=dim(Bi)=dim(Bo)=d, dim(co)=2
% restricted simulation scenario


yalmip('clear')


% C, Cs in I1 O1 I2 O2 F (here F = co) 

Cs = sdpvar(((d*d)^2)*2,((d*d)^2)*2,'symmetric','real');
C = sdpvar(((d*d)^2)*2,((d*d)^2)*2,'symmetric','real');
sdpvar p

F = [];

if isequal(order,[1 2]) % AB
    for i=1:N
        F = F + [linkFAST(kron(J_A(:,:,i),J_B(:,:,i)),Cs,[(d*d)^2 2])==p*linkFAST(kron(J_A(:,:,i),J_B(:,:,i)),SWITCH,[(d*d)^2 2])];
    end
elseif isequal(order,[2 1]) % BA
    for i=1:N
        F = F + [linkFAST(kron(J_B(:,:,i),J_A(:,:,i)),Cs,[(d*d)^2 2])==p*linkFAST(kron(J_A(:,:,i),J_B(:,:,i)),SWITCH,[(d*d)^2 2])];
    end
elseif isequal(order,[1 1]) % AA
    for i=1:N
        F = F + [linkFAST(kron(J_A(:,:,i),J_A(:,:,i)),Cs,[(d*d)^2 2])==p*linkFAST(kron(J_A(:,:,i),J_A(:,:,i)),SWITCH,[(d*d)^2 2])];
    end
elseif isequal(order,[2 2]) % BB
    for i=1:N
        F = F + [linkFAST(kron(J_B(:,:,i),J_B(:,:,i)),Cs,[(d*d)^2 2])==p*linkFAST(kron(J_B(:,:,i),J_B(:,:,i)),SWITCH,[(d*d)^2 2])];
    end    
end

F = F + [Cs>=0,C-Cs>=0,trace(C)==d^2,C==ProjSeqSuperChannel(C,[1 d d d d 2])];


solution = solvesdp(F,-p,sdpsettings('solver','mosek','verbose',1,'cachesolvers',1))


p = double(p);
Cs = double(Cs);
C = double(C);
