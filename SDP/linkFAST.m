%This implements the link product of two matrices A and B, where B is the Choi operator of a linear map from L(C^din) to L(C^dout) and A is an operator in L(C^din). For instance, B may be a quantum channel and A may be a quantum state.
%This code is useful to output C = A*B = \tr_in( (A^T \otimes Id) B) in a FAST way.

%Requires: Nothing
%Input:	A, square matrix din by din
%	B, square matrix din*dout by din*dout
%	d=[din dout]
%Output: C = A*B = \tr_in( (A^T \otimes Id) B) , 

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 01/10/2024

function C = linkFAST(A,B,d)

% A in IN
% B in IN OUT
% C in OUT
% d = [din dout];

din  = d(1);
dout = d(2);

C_vec = kron((conj(A(:)')),eye(dout^2))*PermuteSystems(B(:),[1 3 2 4],[din dout din dout]);

C = reshape(C_vec,[dout dout]);

end
