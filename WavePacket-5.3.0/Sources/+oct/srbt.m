%--------------------------------------------------------------------------
% This function delivers a balancing transform by means of the
% Square Root Balanced Truncation (SRBT) algorithm such that the
% the controllability and observability Gramians become diagonal 
%
% Sigma   = S * W_C * S'
% Sigma   = T' * W_O * T
% Sigma^2 = S * W_C*W_O * T
%
% where T=inv(S). The transformations S,T are found as follows:
%
% (1) Cholesky factorization of Gramians X*X'=WC, and Y*Y'=WO
% for positive definite Gramians performed in the try-block
% Extended to non-positive definite Gramians by setting all negative 
% eigenvalues to zero (done in catch-block)
%
% (2) Singular value decomposition: Y'*X = U*Sigma*V'
%
% (3) Balancing transforms: T=X*V*Sigma^-0.5; S=Sigma^-0.5*U'*Y'
%
% M. S. Tombs, I. Postlethwaite
% Int. J. Control 46(4), pp. 1319-1330, (1987) 
% doi:10.1080/00207178708933971 
% 
% A Laub, M. Heath, C. Paige, R. Ward
% IEEE Trans. Automat. Control 32(2), pp. 115-122, (1987) 
% doi:10.1109/TAC.1987.1104549
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung
%
% see the README file for license details.

function [S, T, Sigma]=srbt(WC, WO)
% threshold = eps('single');

util.disp ('Square Root Balanced Truncation (SRBT)')
util.disp ('======================================')

% Cholesky factorization of W_C
util.disp ('   ')
util.disp ('Trying Cholesky factorization of Gramian W_C = X*X^T')
X = cholesky (WC);

% Cholesky factorization of W_O
util.disp ('   ')
util.disp ('Trying Cholesky factorization of Gramian W_O = Y*Y^T')
Y = cholesky (WO);

% Singular value decomposition should produce a diagonal matrix Sigma 
% of the same dimension as the input matrix, with nonnegative diagonal 
% elements in decreasing order, and unitary matrices U and V.
% Note that the squares of the singular values diag(Sigma) are
% also the eigenvalues of M'*M and M*M'
util.disp ('   ')
util.disp ('Singular value decomposition Y^T*X = U*Sigma*V^T')
M = Y'*X;
[U,Sigma,V] = svd(M);
Sigma = diag(Sigma);
threshold = max(size(M))*eps(norm(M)); % from Matlab's "rank" function
r = nnz(Sigma>=threshold);
util.disp (['Largest singular value                    : ' num2str(max(Sigma))])
util.disp (['Smallest singular value                   : ' num2str(min(Sigma))])
util.disp (['Value of threshold used                   : ' num2str(threshold)])
util.disp (['Dimension of the matrices X and Y         : ' int2str(size(X,1))])
util.disp (['Number of singular values above threshold : ' int2str(r)])
util.disp ( 'Neglecting all singular values below threshold!')
util.disp ('   ')

% Find the balancing transform and its (pseudo-)inverse
% omitting colums/rows corresponding to
% near-zero Hankel singular values 
T = X*V(:,1:r)*diag(Sigma(1:r).^(-1/2));
S = diag(Sigma(1:r).^(-1/2))*(U(:,1:r))'*Y';

end



% This function, loosely speaking, takes the square root of matrix W = X*X'
% If W is symmetric and positive definite, this can be achieved by means of
% Cholesky decomposition in which case X is a lower triangular matrix.
% Alternatively, this "square root" can be obtained as X=V*D^(1/2) where
% D is a diagonal matrix with A's eigenvalues on the main diagonal and
% V is a unitary matrix whose columns are the eigenvectors of A.
% In that case, positive (semi-)definiteness can be enforced by setting
% all eigenvalues below a (zero or small positive) threshold to the value
% of that threshold.
function X = cholesky (W)
threshold = eps('double');
try % for positive definite W
    X=chol(W,'lower');
    util.disp ('Cholesky factorization of Gramian appears to work fine')
catch % for non-positive definite W
    [eigvec,eigval]=eig(W);
    eigval = diag(eigval);
    util.disp ( 'Gramian is not positive definite' )
    util.disp (['Largest eigenvalue                    : ' num2str(max(eigval))])
    util.disp (['Smallest eigenvalue                   : ' num2str(min(eigval))])
    util.disp (['Value of threshold used               : ' num2str(threshold)])
    util.disp (['Number of eigenvalues below threshold : ' int2str(nnz(eigval<threshold))])
    util.disp ( 'Replacing all eigenvalues below threshold!')

    eigval=max(eigval,threshold);
    X=eigvec*sqrt(diag(eigval));
end

end
