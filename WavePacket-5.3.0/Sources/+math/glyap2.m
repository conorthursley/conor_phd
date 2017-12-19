%------------------------------------------------------------------------------
%
% This function determines the solution of the generalized Lyapunov
% equation yielding the matrix X which solves 
%
% A*X + X*A' + sum_d N_d*X*N_d' + B*B' = 0,  for tflag==false  
% A'*X + X*A + sum_d N_d'*X*N_d + B'*B = 0,  for tflag==true
%
% using the biconjugate gradient method (with preconditioner)
% as suggested by Tobias Breiten from TU Graz, Austria
% 
% One can consider the generalized Lyapunov equation as a system of linear
% equations A~*x=b where x is a vectorization of the unknown matrix X and 
% the action of A~ on x corresponds to the first two terms in the equations
% given above. This is routinely solved by the biconjugate gradient method.
% Using a preconditioner M one effectively solves the modified system 
% inv(M)*A~*x = inv(M)*b. If one chooses M such that M*x=b is "close" to the 
% original system, then inv(M)*A~ will be close to the identity matrix thus
% rendering the modified system (hopefully!) well-conditioned.
% Here M is derived from solution of the standard Lyapunov equation.
%
% params.conv_tol  = Convergence tolerance: default 1e-9
% params.max_iter  = Maximal number of iterations; default 50
%
% The  above generalized Lyapunov equations arise, e.g., in the context of 
% bilinear control systems, where they are used to calculate the Gramians.
% See, for example
% 
% L. Zhang, J. Lam: Automatica 38, 205 (1988)
% doi:10.1016/S0005-1098(01)00204-7
%
% Z. Bai, D. Skoogh, Lin. Alg. Appl. 415, 406 (2006)
% doi:10.1016/j.laa.2005.04.032
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014 Tobias Breiten
%               2015 Burkhard Schmidt
%
% see the README file for license details.

function X = glyap2(A,B,N,tflag,params) 
 
% Check/set parameters
n=size(A,1);
util.disp ('Solving generalized Lyapunov equation by biconjugate gradients')
util.disp (['for dimensionality = ' int2str(n)])

if ~isfield(params,'conv_tol')
    params.conv_tol = 1e-10;
end
util.disp (['Convergence tolerance = ' num2str(params.conv_tol)])

if ~isfield(params,'max_iter')
    params.max_iter = 10;
end
util.disp (['Max. number of iterations = ' int2str(params.max_iter)])

% Optionally transposing all input vectors/matrices
if tflag
    A = A';
    for d=1:length(N) 
        N{d} = N{d}';
    end
    for d=1:length(B) 
        B{d} = B{d}';
    end
end

% Construct right-hand-side vector b
b = zeros(n);
for d=1:length(B)
    b = b - B{d}*B{d}';
end

% Call  bilinear conjugate gradient method to solve A~*x=b
X = bicg( ...
     @(x,transp_flag)glyap_afun(x,n,A,N,transp_flag), ...
     reshape(b,n^2,1), ...
     params.conv_tol, ...
     params.max_iter, ...
     @(x,transp_flag)glyap_mfun(x,n,A,transp_flag));
X = reshape(X,n,n); 

% Finally check (error norm of) deviation 
Deviate = A*X + X*A' - b;
for d=1:length(N)
    Deviate = Deviate+N{d}*X*N{d}';
end
util.disp(['GLYAP2 error norm = ' num2str(norm(Deviate))]);
util.disp ('   ')
end

% A~-function for bilinear conjugate gradient method to solve A~*x=b
% afun(x,'notransp') returns A~*x and afun(x,'transp') returns A~'*x
function y = glyap_afun(x,n,A,N,transp_flag)

X = reshape(x,n,n);                      % Reshape: vector ==> matrix

if strcmp(transp_flag,'transp')          % Y = A'*X + X*A + sum (N'*X*N)
    Y = A'*X + X*A;
    for d = 1:length(N)
        Y = Y + N{d}'*X*N{d};
    end    
elseif strcmp(transp_flag,'notransp')    % Y = A*X + X*A' + sum (N*X*N')
    Y = A*X + X*A';
    for d = 1:length(N)
        Y = Y + N{d}*X*N{d}';
    end
end

y = Y(:);                            % Reshape: matrix ==> vector
end

% M-function for pre-conditioning the bilinear conjugate gradient method
% mfun(x,'notransp') returns M\x and mfun(x,'transp') returns M'\x.
% here: using solution of the standard Lyapunov equation
function y = glyap_mfun(x,n,A,transp_flag)

X = reshape(x,n,n);                     % Reshape: vector ==> matrix

if strcmp(transp_flag,'transp')
    Y = lyap(A',X);
elseif strcmp(transp_flag,'notransp')
    Y = lyap(A,X);
end

y = Y(:);                                % Reshape: matrix ==> vector
end
