%------------------------------------------------------------------------------
%
% This function determines the solution of the generalized Sylvester
% equation yielding the matrix X which solves 
%
% A*X + X*Ar' + sum_d N_d*X*Nr_d' + B*Br' = 0,  for tflag==false  
% A'*X + X*Ar + sum_d N_d'*X*Nr_d - B'*Br = 0,  for tflag==true
%
% Note the minus sign in the second equation!
%
% using the biconjugate gradient method (with preconditioner)
% as suggested by Tobias Breiten from TU Graz, Austria
%
% One can consider the generalized Sylvester equation as a system of linear
% equations A~*x=b where x is a vectorization of the unknown matrix X and 
% the action of A~ on x corresponds to the first two terms in the equations
% given above. This is routinely solved by the biconjugate gradient method.
% Using a preconditioner M one effectively solves the modified system 
% inv(M)*A~*x = inv(M)*b. If one chooses M such that M*x=b is "close" to the 
% original system, then inv(M)*A~ will be close to the identity matrix thus
% rendering the modified system (hopefully!) well-conditioned.
% Here M is derived from solution of the standard Sylvester equation.
%
% params.conv_tol  = Convergence tolerance: default 1e-9
% params.max_iter  = Maximal number of iterations; default 50
%
% Interpolation-Based H2-Model Reduction of Bilinear Control Systems;
% Benner, Peter; Breiten, Tobias;
% SIAM Journal on Matrix Analysis and Applications : 
% Vol. 33, No. 3, 859-885; 2012.     doi:10.1137/110836742
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014 Tobias Breiten
%               2015 Burkhard Schmidt
%
% see the README file for license details.

function X = gsylv2 (A,B,N,Ar,Br,Nr,tflag,params)
n=size(A, 1);
r=size(Ar,1);

if ~isfield(params,'conv_tol')
    params.conv_tol = 1e-10;
end
% util.disp (['Convergence tolerance = ' num2str(params.conv_tol)])

if ~isfield(params,'max_iter')
    params.max_iter = 10;
end
% util.disp (['Max. number of iterations = ' int2str(params.max_iter)])

% Optionally transposing all input vectors/matrices
if tflag
    A  = A';
    Ar = Ar';
    for d=1:length(N) 
        N {d} = N {d}';
        Nr{d} = Nr{d}';
    end
    for d=1:length(B) 
        B {d} = B {d}';
        Br{d} = Br{d}';
    end
end

% Construct right-hand-side vector b
b = zeros(n,r);
for d=1:length(B)
    b = b - B{d}*Br{d}';
end
if tflag
    b = - b;
end

% Call  bilinear conjugate gradient method to solve A~*x=b
X = bicg( ...
    @(x,transp_flag)gsylv_afun(x,n,r,A,N,Ar,Nr,transp_flag), ...
    reshape(b,n*r,1), ...
    params.conv_tol, ...
    params.max_iter, ...
    @(x,transp_flag)gsylv_mfun(x,n,r,A,Ar,transp_flag));
X = reshape(X,n,r);

end

% A~-function for bilinear conjugate gradient method to solve A~*x=b
% afun(x,'notransp') returns A~*x and afun(x,'transp') returns A~'*x
function y = gsylv_afun(x,n,r,A,N,Ar,Nr,transp_flag)

X = reshape(x,n,r);                      % Reshape: vector ==> matrix

if strcmp(transp_flag,'transp')          % y = A'*x
    Y = A'*X + X*Ar;
    for d = 1:length(N)
        Y = Y + N{d}'*X*Nr{d};
    end
elseif strcmp(transp_flag,'notransp')    % y = A*x
    Y = A*X + X*Ar';
    for d = 1:length(N)
        Y = Y + N{d}*X*Nr{d}';
    end
end

y = Y(:);                                % Reshape: matrix ==> vector
end

% M-function for pre-conditioning the bilinear conjugate gradient method
% mfun(x,'notransp') returns M\x and mfun(x,'transp') returns M'\x.
% here: using solution of the standard Sylvester equation
function y = gsylv_mfun(x,n,r,A,Ar,transp_flag)

X = reshape(x,n,r);                     % Reshape: vector ==> matrix

if strcmp(transp_flag,'transp')
    Y = math.sylv(A',Ar,X);
elseif strcmp(transp_flag,'notransp')
    Y = math.sylv(A,Ar',X);
end

y = Y(:);                                % Reshape: matrix ==> vector
end