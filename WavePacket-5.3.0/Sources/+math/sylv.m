function [X] = sylv (A,H,M,Trans)

%--------------------------------------------------------------------------
% Solving the Sylvester equation
%                AX+XH+M=0 
% with H small and dense and A large and sparse. 
% If Trans=="T" the transposed equation 
%               A'X+XH'+M=0
% is solved. If Trans!="T" or it is not set, the 
% normal equation is solved.
%
% Observation: For the cases where one of the matrices is much smaller than
% the other, this code is (much!) faster than function 'lyap' contained in
% Matlab's Control System Toolbox
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C)  2011 Martin Koehler MPI Magdeburg

if ( nargin < 4)
    Trans = 'N';
end

q = size(H,1);
if isreal(H)
    isr = 1;
else
    isr = 0;
end

[U,S]=schur(H,'complex');
Mtilde = M*U;
Xtilde = zeros(size(A,1),q);
I=speye(size(A));

if (Trans=='N')
    for j = 1:q
        rhs=-Mtilde(:,j);
        for i=1:j-1
            rhs = rhs - S(i,j)*Xtilde(:,i);
        end
        Xtilde(:,j)=(A+S(j,j)*I)\rhs;
    end
else
    %S=S';
    for jj=1:q
        j=q-jj+1;
        rhs=-Mtilde(:,j);
         for i=j+1:q
            rhs = rhs - conj(S(j,i))*Xtilde(:,i);
        end
        Xtilde(:,j)=(A+S(j,j)*I)'\rhs;
    end
end
            
X=Xtilde*U';

if isr 
    X=real(X);
end
eigH=diag(S);

end

