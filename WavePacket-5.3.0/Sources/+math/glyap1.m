%------------------------------------------------------------------------------
%
% This function determines the solution of the generalized Lyapunov
% equation yielding the matrix X which solves 
%
% A*X + X*A' + sum_d N_d*X*N_d' + B*B' = 0,  for tflag==false  
% A'*X + X*A + sum_d N_d'*X*N_d + B'*B = 0,  for tflag==true
%
% using an iterative scheme introduced inde 
% E. Wachspress, Appl. Math. Lett. 1, 87 (1988)
% 
% params.dev_max  = Initial maximal deviation: default 1e10
% params.conv_tol = Convergence tolerance: default 1e-10
% params.max_iter = Maximal number of iterations; default 100
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
% Copyright (C) 2011 Boris Schaefer-Bung
%               2014 Burkhard Schmidt
%
% see the README file for license details.

function X = glyap1(A,B,N,tflag,params) 
 
% Check/set parameters
n=size(A,1);
util.disp ('Solving generalized Lyapunov equation by iteration method')
util.disp (['for dimensionality = ' int2str(n)])

if ~isfield(params,'dev_max')
    params.dev_max = 1e10;
end
util.disp (['Initial maximal deviation = ' num2str(params.dev_max)])

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

% Construct vector B*B'
BB = zeros(n);
for d=1:length(B)
    BB = BB + B{d}*B{d}';
end

% Initial guess: Start the iteration from solution of the linear problem
% i.e., the solution of the ordinary Lyapunov equation (without matrices N)
X = lyap(A,BB);

% Iteration scheme for generalized Lyapunov equation
util.disp('Iteration sum|W_i-W_{i-1}| max|W_i-W_{i-1}| sum|W_i|    rank(W_i)')
for iter=1:params.max_iter
    Xold=X;
    NXN = zeros(n);
    for d=1:length(N)
        NXN = NXN + N{d}*Xold*N{d}';
    end
    
    % Standard Lyapunov equation solver from Matlab's control tool box
    X=lyap(A,NXN+BB);
    X=0.5*(X+X'); % Enforce symmetry
    
    % Calculate errors and check for convergence
    dev_sum=sum(abs(X(:)-Xold(:)));
    dev_max=max(abs(X(:)-Xold(:)));
    sum_mat=sum(abs(X(:)));
    rkcrhs=rank(X); % output rank
    util.disp(['    ' num2str(iter, '%3.0d') '      ' num2str(dev_sum, '%11.4e')...
        '       ' num2str(dev_max, '%11.4e') '       ' ...
            num2str(sum_mat, '%11.4e') '  ' num2str(rkcrhs)])
    if ((dev_sum < params.conv_tol) && (dev_sum/sum_mat < params.conv_tol))
        util.disp(['Convergence reached after ' num2str(iter) ' iterations'] )
        break
    end
    if (params.dev_max < dev_max)
        util.disp(['Smallest deviation after ' num2str(iter-1) ' iterations'] )
        X=Xold;
        break
    end
    params.dev_max=dev_max;
end

% Finally check (error norm of) deviation 
Deviate = A*X + X*A' + BB;
for d=1:length(N)
    Deviate = Deviate + N{d}*X*N{d}';
end
util.disp(['GLYAP1 terminated after ' int2str(iter) ' iterations with error norm = ' num2str(norm(Deviate))]);
util.disp ('   ')
end

