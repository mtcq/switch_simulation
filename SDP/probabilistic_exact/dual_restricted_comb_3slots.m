% Author: Jessica Bavaresco, jessicabavaresco.github.io
% Requires: QETLAB and all linkFAST.m
% Last update: 01 Oct 2024

function [G,Ri,p] = dual_restricted_comb_3slots(J_A,J_B,d,order)

% Compared to the paper, here we use a relabelling (i,j) -> i of variables Rij -> Ri and {J^A_i, J^B_j} -> {J^A_i, J^B_j}. 

% INPUT: J_A and J_B are N pairs of input channels, with size(J_A) = size(J_B) = [d^2 d^2 N]
%        d the dimension of the target systems 

% OUTPUT: G (called \Gamma in the paper) and Ri are the optimal solutions of the dual problem, with size(Ri) = [2 2 N]
%         p is the maximum probability of success  

N = size(J_A,3); 

for i=1:N
    J_A(:,:,i) = transpose(J_A(:,:,i));
    J_B(:,:,i) = transpose(J_B(:,:,i));
end

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


% G in I1 O1 I2 O2 I3 O3 F (here F = co)
% Ri in F

G  = sdpvar(((d*d)^3)*2,((d*d)^3)*2,'hermitian','complex');
Ri = sdpvar(2,2,N,'hermitian','complex');

t=0;
RD = zeros(((d*d)^3)*2,((d*d)^3)*2);

if isequal(order,[1 2 1]) % ABA
    for i=1:N
        t = t + trace(Ri(:,:,i)*PartialTrace(SWITCH*Tensor(J_A(:,:,i),J_B(:,:,i),eye(2)),1,[(d*d)^2 2])); 
        RD  = RD + Tensor(J_A(:,:,i),J_B(:,:,i),J_A(:,:,i),Ri(:,:,i));
    end
elseif isequal(order,[2 1 2]) % BAB
    for i=1:N
        t = t + trace(Ri(:,:,i)*PartialTrace(SWITCH*Tensor(J_A(:,:,i),J_B(:,:,i),eye(2)),1,[(d*d)^2 2])); 
        RD  = RD + Tensor(J_B(:,:,i),J_A(:,:,i),J_B(:,:,i),Ri(:,:,i));
    end
elseif isequal(order,[1 1 2]) % AAB
    for i=1:N
        t = t + trace(Ri(:,:,i)*PartialTrace(SWITCH*Tensor(J_A(:,:,i),J_B(:,:,i),eye(2)),1,[(d*d)^2 2])); 
        RD  = RD + Tensor(J_A(:,:,i),J_A(:,:,i),J_B(:,:,i),Ri(:,:,i));
    end
elseif isequal(order,[2 2 1]) % BBA
    for i=1:N
        t = t + trace(Ri(:,:,i)*PartialTrace(SWITCH*Tensor(J_A(:,:,i),J_B(:,:,i),eye(2)),1,[(d*d)^2 2])); 
        RD  = RD + Tensor(J_B(:,:,i),J_B(:,:,i),J_A(:,:,i),Ri(:,:,i));
    end
elseif isequal(order,[2 1 1]) % BAA
    for i=1:N
        t = t + trace(Ri(:,:,i)*PartialTrace(SWITCH*Tensor(J_A(:,:,i),J_B(:,:,i),eye(2)),1,[(d*d)^2 2])); 
        RD  = RD + Tensor(J_B(:,:,i),J_A(:,:,i),J_A(:,:,i),Ri(:,:,i));
    end
elseif isequal(order,[1 2 2]) % ABB
    for i=1:N
        t = t + trace(Ri(:,:,i)*PartialTrace(SWITCH*Tensor(J_A(:,:,i),J_B(:,:,i),eye(2)),1,[(d*d)^2 2])); 
        RD  = RD + Tensor(J_A(:,:,i),J_B(:,:,i),J_B(:,:,i),Ri(:,:,i));
    end
elseif isequal(order,[1 1 1]) % AAA
    for i=1:N 
        t = t + trace(Ri(:,:,i)*PartialTrace(SWITCH*Tensor(J_A(:,:,i),J_A(:,:,i),eye(2)),1,[(d*d)^2 2])); 
        RD  = RD + Tensor(J_A(:,:,i),J_A(:,:,i),J_A(:,:,i),Ri(:,:,i));
    end
elseif isequal(order,[2 2 2]) % BBB
    for i=1:N
        t = t + trace(Ri(:,:,i)*PartialTrace(SWITCH*Tensor(J_B(:,:,i),J_B(:,:,i),eye(2)),1,[(d*d)^2 2])); 
        RD  = RD + Tensor(J_B(:,:,i),J_B(:,:,i),J_B(:,:,i),Ri(:,:,i));
    end    
end

F = [t==1,G-RD>=0,G>=0,G==ProjSeqSuperChannel(G, [d d d d d d 2 1])];

p = (1/((d^3)*2))*trace(G);


solution = solvesdp(F,p,sdpsettings('solver','mosek','verbose',1,'cachesolvers',1));


G  = double(G);
Ri = double(Ri);
p  = double(p);



