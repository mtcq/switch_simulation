%Author: Marco T??lio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Requires: nothing
%Last update: 03/03/2024

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
