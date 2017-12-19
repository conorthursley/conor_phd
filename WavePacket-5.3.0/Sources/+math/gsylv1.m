%--------------------------------------------------------------------------
% 
% This function determines the solution of the generalized Sylvester
% equation yielding the matrix X which solves 
%
% A*X + X*Ar' + sum_d N_d*X*Nr_d' + B*Br' = 0,  for tflag==false  
% A'*X + X*Ar + sum_d N_d'*X*Nr_d - B'*Br = 0,  for tflag==true
%
% Note the minus sign in the second equation!
%
% using the iterative method (with preconditioner)
% as suggested by Tobias Breiten from TU Graz, Austria
% 
% params.conv_tol = Convergence tolerance: default 1e-8
% params.max_iter = Maximal number of iterations; default 10
%
% Typically, this converges after 3-4 steps
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

function X = gsylv1 (A,B,N,Ar,Br,Nr,tflag,params)
n = size(A, 1);
r = size(Ar,1);

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

% Construct vector B*Br'
BBr = zeros(n,r);
for d=1:length(B)
    BBr = BBr + B{d}*Br{d}';
end
if tflag
    BBr = - BBr;
end

% Initial guess: Start the iteration from solution of the linear problem
% i.e., the solution of the ordinary Sylvester equation (without matrices N)
X = math.sylv(A,Ar',BBr);
% X = sylvester(full(A),Ar',-BBr);
    
% Iteration scheme for generalized Sylvester equation
err = 1;
iter = 0;
while(err > params.conv_tol && iter < params.max_iter)
    Xold = X;
    NXNr = zeros(n,r);
    for d = 1:length(N)
        NXNr = NXNr + N{d}*Xold*Nr{d}';
    end
    
    % Magdeburg solver is (much!) faster than Matlab's function 'lyap'
    X = math.sylv(A,Ar',NXNr+BBr);
    
    % Calculate error and check for convergence
    err = norm(X-Xold,'fro')/norm(Xold,'fro');
    % util.disp ([int2str(iter) ': error = ' num2str(e)])
    iter = iter + 1;
end

% Finally check (error norm of) deviation 
Deviate = A*X + X*Ar' + BBr;
for d=1:length(N)
    Deviate = Deviate + N{d}*X*Nr{d}';
end
util.disp(['GSYLV1 terminated after ' int2str(iter) ' iterations with error norm = ' num2str(norm(Deviate))]);
end
