%--------------------------------------------------------------------------
%
% H2 optimal model reduction
%
% Solving a generalized Sylvester equation by BIRKA : The Bilinear
% Iterative Rational Krylov Algorithm (BIRKA) aims at the construction 
% of bilinear reduced-order models that are (locally) optimal with respect 
% to the bilinear H2-norm. Based on the idea of iterative correction, 
% the algorithm can be seen as a suitable extension of its linear analogue 
% (IRKA). Locally H2-optimal reduced-order models fulfill an abstract 
% Hermite-type interpolation condition for the Volterra series of the 
% bilinear control system. Similar as in the linear case, these optimality 
% conditions can also be stated in terms of (generalized) linear matrix 
% equations. This in turn allows to construct the required projection 
% subspaces as solutions to certain generalized Sylvester equations. 
% 
%
% A, B, N, C matrices and initial/equilibrium state vectors are read from
% "filename.mat" where filename typically is 'lvne' or 'tdse'. In the end, 
% the reduced matrices and vectors are written to "filename_r<n>.mat". 
% Input parameter 'dim' specifies the dimension of the truncated system.
%
% reduce.H2model.max_iter:
% Since BIRKA is an iterative procedure, we need some iteration cycle 
% that should not exceed the number max_iter. For non-converging 
% initialization, max_iter prevents from an infinite iteration cycle. 
% Typical values are max_iter = 50,100,200 (depending on the size of 
% the original model).
%
% reduce.H2model.conv_tol:
% In numerical computations, the iteration cycle will never yield a 
% convergence tolerance which is exactly 0. We therefore have to specify 
% a value conv_tol which can be seen as a numerical perturbation that we  
% expect not to influence the quality of the computed model significantly.
%
% Interpolation-Based H2-Model Reduction of Bilinear Control Systems;
% Benner, Peter; Breiten, Tobias;
% SIAM Journal on Matrix Analysis and Applications : 
% Vol. 33, No. 3, 859–885; 2012.     doi:10.1137/110836742
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2013    Tobias Breiten and Peter Benner, MPI Magdeburg 
%               2014-16 Burkhard Schmidt, FU Berlin

function qm_H2model(filename, dim)

global bilinear reduce 

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

util.disp (' ')
util.disp ('--------------------------------------------------------------------------------')
util.disp (' Interpolation based H2 model reduction (c) 2013, Breiten/Benner, MPI Magdeburg ')
util.disp (' Reducing A, B, N, C matrices; also initial/equilibrium state vectors, x_i, x_e ')
util.disp (' http://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_H2model   ')
util.disp ('--------------------------------------------------------------------------------')
util.disp (' ')
util.disp(['reducing to ' int2str(dim) ' modes:'])
util.disp (' ')

% Set all parameters to default values if not specified by the user 
if ~isfield (reduce,'H2model')
    reduce.H2model = [];
end

if ~isfield (reduce,'gsylv')
    reduce.gsylv=[];
end

if ~isfield (reduce.H2model,'conv_tol')
    reduce.H2model.conv_tol=1e-6;
end
util.disp(['Convergence tolerance for BIRKA iteration: ' num2str(reduce.H2model.conv_tol)])

if ~isfield (reduce.H2model,'max_iter')
    reduce.H2model.max_iter=100;
end
util.disp(['Maximal number of BIRKA iterations: ' int2str(reduce.H2model.max_iter)])

if ~isfield (reduce.H2model,'method')
    reduce.H2model.method='iter';
end
util.disp(['Method for solving generalized Sylvester equation: ' reduce.H2model.method])

if ~isfield(reduce.H2model,'BN_scale')
    reduce.H2model.BN_scale = 1;
end
util.disp(['Scaling B-vector, N-matrix: ' num2str(reduce.H2model.BN_scale)])

% Load full matrices from data file  
% downscale B and N such that generalized Lyapunov equations are solvable
util.disp (['Loading bilinear system matrices from file: ' filename '.mat'])
load(filename)
n = size(bilinear.A,1);
util.disp (['Dimensionality of original system: ' int2str(n)])

A = bilinear.A; 

B = cell(length(bilinear.B),1);
for d=1:length(bilinear.B)
    B{d} = bilinear.B{d} / reduce.H2model.BN_scale;
end

N = cell(length(bilinear.N),1);
for d=1:length(bilinear.N)
    N{d} = bilinear.N{d} / reduce.H2model.BN_scale;
end

C = cell(length(bilinear.C),1);
for d=1:length(bilinear.C)
    C{d} = bilinear.C{d};
end
util.disp ('   ')

% Check stability: real part of all eigenvalues should be negative
oct.check_stable ('original system matrix A', A)
util.disp (' ')

%% Two choices of stabilizing A matrix
util.disp (' ')
util.disp ('Stabilizing A matrix')
if ~isfield(reduce.H2model,'A_stable')
    reduce.H2model.A_stable = 'ssu';
end
switch lower(reduce.H2model.A_stable)
    case 'ssu' 
        % splitting stable from unstable part i.e. where Re(eig(A)) >= 0
        % consider only the orthogonal complement of zero eigenspace
        % generalization: split all eigenvalues with real part below
        % specified threshold

        util.disp ('Method ''ssu'' : splitting stable from unstable part')
        if ~isfield (reduce.H2model,'A_split')
            reduce.H2model.A_split = 1;
        end
        n1 = reduce.H2model.A_split; 
        util.disp (['Dimension of unstable part : ',int2str(n1)])
        util.disp (' ')
        
        % Calculate all eigenvalues and eigenvectors (columns of U) of A.
        % In general (T>0, with decoherence), A is complex & non-symmetric
        [U,D] = eig(full(A));
        [~, ind] = sort(real(diag(D)), 'descend');
        U = U(:,ind);
                         
        % Transformation and splitting of A, B, N, C into eigenvector basis
        A  = U\A*U;             % transformation; should be diagonal matrix
        A  = util.sparse(A);    % reinforce sparsity of A
        Au = A(n1+1:n,   1:n1); % unstable part
        A  = A(n1+1:n,n1+1:n ); % stable part
        
        for d=1:length(B)
            B{d} = U\B{d};       % transformation
            B{d} = B{d}(n1+1:n); % stable part
        end
         
        Nu = cell(length(N),1);
        for d=1:length(N) 
            N{d}  = U\N{d}*U;             % transformation
            N{d}  = util.sparse(N{d});    % reinforce sparsity of A
            Nu{d} = N{d}(n1+1:n,   1:n1); % unstable part
            N{d}  = N{d}(n1+1:n,n1+1:n ); % stable part
        end
        
        for d=1:length(C)
            C{d} = C{d}*U;       % transformation
            C{d} = C{d}(n1+1:n); % stable part 
        end
        
        % Rho coupling terms act as additional control field
        % Does it really make things any better?
        if ~isfield (reduce.H2model,'acf_couple')
           reduce.H2model.acf_couple=0; 
        end
        if reduce.H2model.acf_couple
            nBg = length(B);
            B{nBg+1} = Au;
            for d=1:length(N) 
                B{nBg+1+d} = Nu{d};
            end
        end
        
    case 'evs' 
        util.disp('Method ''evs'': eigenvalue shift of A')
        
        % value of shift
        if ~isfield(reduce.H2model,'A_shift')
            reduce.H2model.A_shift = 1e-6;
        end
        util.disp (['Value of shift : ' num2str(reduce.H2model.A_shift)])
        util.disp (' ')
        
        % shift diagonal values of A-matrix
        A = A - reduce.H2model.A_shift * speye(n);
                
    otherwise
        util.error (['Wrong choice of stabilization method : ' reduce.H2model.A_stable])
        
end

% Check stability of the splitted (SSU) or shifted (EVS) system
switch lower(reduce.H2model.A_stable)
    case 'ssu' 
        oct.check_stable ('splitted system matrix A', A)
    case 'evs' 
        oct.check_stable ('shifted system matrix A', A)
end
util.disp (' ')

% Initialization parameters for the iteration cycle. 
% As long as a reasonable reduced-order model is not known, 
% one might just as well use a random initialization
rng('default'); % set random number generator as if you restarted MATLAB
                   Ar    = rand(dim);       % full matrix
for d=1:length(B); Br{d} = rand(dim,1); end % column vector(s)
for d=1:length(N); Nr{d} = rand(dim);   end % full matrix(s)
for d=1:length(C); Cr{d} = rand(1,dim); end % row vector(s)

%% Start BIRKA iteration loop
Spec = eig(Ar);
Spec = -sort(Spec);
err = inf;
iter = 0;
while(err > reduce.H2model.conv_tol && iter < reduce.H2model.max_iter)
    Spec_old = Spec;
    
    % Solve generalized Sylvester equations
    switch lower (reduce.H2model.method)
        case 'iter'
            X = math.gsylv1(A,B,N,Ar,Br,Nr,false,reduce.gsylv);
            Y = math.gsylv1(A,C,N,Ar,Cr,Nr,true ,reduce.gsylv);
        case 'bicg'
            X = math.gsylv2(A,B,N,Ar,Br,Nr,false,reduce.gsylv);
            Y = math.gsylv2(A,C,N,Ar,Cr,Nr,true ,reduce.gsylv);
        otherwise
            util.error (['Wrong choice of GSYLV solver method : ' reduce.H2model.method])    
    end
    
    % Orthogonal-triangular decomposition A = Q * R where A is m-by-n
    % where Q is m-by-m (unitary matrix Q)
    % where R is m-by-n (upper triangular matrix)
    [T,~] = qr(X,0);
    [W,~] = qr(Y,0);
    
    % Refining reduced matrices
    S = (W'*T)\W';
                       Ar    = S*A   *T;
    for d=1:length(B); Br{d} = S*B{d}  ; end
    for d=1:length(N); Nr{d} = S*N{d}*T; end
    for d=1:length(C); Cr{d} =   C{d}*T; end
    
    % Spectrum of A
    Spec = eig(Ar);
    Spec = -sort(Spec);
    iter = iter + 1;
    
    % Check change in spectrum of A for convergence
    err = norm(Spec-Spec_old)/norm(Spec_old);
    util.disp([int2str(iter) '-th BIRKA step: conv. crit. = ' num2str(err)])
    
end
redmat.S=S;
redmat.T=T;

%% Add extra columns, rows to S, V matrices to compensate for the SSU splitting
if strcmpi(reduce.H2model.A_stable,'ssu')
    dimt=size(redmat.T);
    redmat.S=[eye(n1) zeros(n1,dimt(1)); zeros(dimt(2),n1) redmat.S];
    redmat.S=redmat.S/U;
    redmat.T=[eye(n1) zeros(n1,dimt(2)); zeros(dimt(1),n1) redmat.T];
    redmat.T=U*redmat.T; 
end

%% Constructing reduced matrices
% Check whether T is the (pseudo-)inverse of S 
util.pseudo(redmat.S,redmat.T)

% Transform A,B,N,C
redmat.A = redmat.S*bilinear.A*redmat.T;
for d=1:length(bilinear.B)
    redmat.B{d} = redmat.S*bilinear.B{d}; 
end
for d=1:length(bilinear.N)
    redmat.N{d} = redmat.S*bilinear.N{d}*redmat.T;
end
for d=1:length(bilinear.C)
    redmat.C{d} = bilinear.C{d}*redmat.T;
    redmat.Q{d} = bilinear.Q{d};
end

% State vectors (x, transformed) and output (y, unchanged)
redmat.x.initial = redmat.S*bilinear.x.initial;
redmat.x.equilib = redmat.S*bilinear.x.equilib;
redmat.y.initial=            bilinear.y.initial;
redmat.y.equilib=            bilinear.y.equilib;

% Labels of control targets and title
redmat.label=bilinear.label;  % labels of control targets
redmat.title=[bilinear.title ', H2 error reduction to ']; % simulation title

% Save reduced matrices to data file
clear global bilinear
global bilinear
bilinear = redmat;
util.disp(['Saving reduced matrices A, B, N, C and vectors x,y to file : ' filename '_h' int2str(dim) '.mat'])
save([filename '_h' int2str(dim)],'bilinear');

% Plot spectrum of A 
oct.spec_A (mfilename)

% Output clock/date/time
util.clock;

end

