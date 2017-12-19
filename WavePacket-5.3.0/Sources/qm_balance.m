%------------------------------------------------------------------------------
%
% Balancing transformation of controllability and observability Gramians:
% A, B, N, C matrices and initial/equilibrium state vectors are read from
% "filename.mat" where filename typically is 'lvne' or 'tdse'. In the end,
% the balanced matrices and vectors are written to file "filename_b.mat".
%
% For more details see
% B. Schaefer-Bung, C. Hartmann, B. Schmidt, and Ch. Schuette
% J. Chem. Phys. 135, 014112 (2011)
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung
%               2013-16 Burkhard Schmidt

function qm_balance(filename)

global balanced bilinear reduce

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

util.disp (' ')
util.disp ('------------------------------------------------------------------------------')
util.disp (' Balancing transformation of A, B, N, C matrices; also x_i, x_e state vectors ')
util.disp (' http://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_balance')
util.disp ('------------------------------------------------------------------------------')
util.disp (' ')

% Set all parameters to default values if not specified by the user 
if ~isfield (reduce,'balance')
    reduce.balance=[];
end

if ~isfield (reduce,'glyap')
    reduce.glyap=[];
end

if ~isfield (reduce.balance,'method')
    reduce.balance.method='iter';
end
util.disp(['Method for solving generalized Lyapunov equation: ' reduce.balance.method])

if ~isfield (reduce.balance,'transform')
    reduce.balance.transform = 'srbt';
end
util.disp(['Method of balancing transformation: ' reduce.balance.transform])

if ~isfield(reduce.balance,'BN_scale')
    reduce.balance.BN_scale = 1;
end
util.disp(['Scaling B-vector, N-matrix: ' num2str(reduce.balance.BN_scale)])

% Load full matrices from data file  
% downscale B and N such that generalized Lyapunov equations are solvable
util.disp (['Loading bilinear system matrices from file: ' filename '.mat'])
load(filename)
n = size(bilinear.A,1);
util.disp (['Dimensionality of original system: ' int2str(n)])

A = bilinear.A;

B = cell(length(bilinear.B),1);
for d=1:length(bilinear.B)
    B{d} = bilinear.B{d}/reduce.balance.BN_scale;
end

N = cell(length(bilinear.N),1);
for d=1:length(bilinear.N)
    N{d} = bilinear.N{d}/reduce.balance.BN_scale;
end

C = cell(length(bilinear.C),1);
for d=1:length(bilinear.C)
    C{d}=bilinear.C{d};
end
util.disp ('   ')

% Check stability: real part of all eigenvalues should be negative
oct.check_stable ('original system matrix A', A)
% oct.check_stable ('original system matrix A+N*N^T', A, N)
util.disp (' ')

%% Two choices of stabilizing A matrix
util.disp ('Stabilizing A matrix')
if ~isfield(reduce.balance,'A_stable')
    reduce.balance.A_stable = 'ssu';
end
switch lower(reduce.balance.A_stable)
    case 'ssu' 
        % splitting stable from unstable part i.e. where Re(eig(A)) >= 0
        % consider only the orthogonal complement of zero eigenspace
        % generalization: split all eigenvalues with real part below
        % specified threshold
        
        util.disp ('Method ''ssu'' : splitting stable from unstable part')     
        if ~isfield (reduce.balance,'A_split')
            reduce.balance.A_split = 1;
        end
        n1 = reduce.balance.A_split; 
        util.disp (['Dimension of unstable part : ',int2str(n1)])
        util.disp (' ')
        
        % Calculate all eigenvalues and eigenvectors (columns of U) of A.
        % Note that in general (T>0, LvNE), A is complex and non-symmetric
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
        if ~isfield (reduce.balance,'acf_couple')
           reduce.balance.acf_couple=0; 
        end
        if reduce.balance.acf_couple
            nBg = length(B);
            B{nBg+1} = Au;
            for d=1:length(N) 
                B{nBg+1+d} = Nu{d};
            end
        end
        
    case 'evs' 
        util.disp('Method ''evs'': eigenvalue shift of A')
        
        % value of shift
        if ~isfield(reduce.balance,'A_shift')
            reduce.balance.A_shift = 1e-6;
        end
        util.disp (['Value of shift : ' num2str(reduce.balance.A_shift)])
        util.disp (' ')
        
        % shift diagonal values of A-matrix
        A = A - reduce.balance.A_shift * speye(n);
                
    otherwise
        util.error (['Wrong choice of stabilization method : ' reduce.balance.A_stable])
        
end

% Check stability of the splitted (SSU) or shifted (EVS) system
switch lower(reduce.balance.A_stable)
    case 'ssu' 
        oct.check_stable ('splitted system matrix A', A)
        % oct.check_stable ('splitted system matrix A+N*N^T', A, N)
    case 'evs' 
        oct.check_stable ('shifted system matrix A', A)
        % oct.check_stable ('shifted system matrix A+N*N^T', A, N)
end
util.disp (' ')
  
%% Controllability Gramian from generalized Lyapunov equation
util.disp('------------------------------------');
util.disp('Solve generalized Lyapunov equation ');
util.disp('to get controllability Gramian W_C  ');
util.disp('A*WC + WC*Ap + N*WC*Np + B*Bp = 0   ');
switch lower (reduce.balance.method)
    case 'iter'
        WC = math.glyap1(A,B,N,false,reduce.glyap);
    case 'bicg'
        WC = math.glyap2(A,B,N,false,reduce.glyap);
    otherwise
        util.error (['Wrong choice of GLYAP solver method : ' reduce.balance.method])
end

%% Observability Gramian from generalized Lyapunov equation
util.disp('------------------------------------');
util.disp('Solve generalized Lyapunov equation ');
util.disp('to get observability Gramian W_O    ');
util.disp('Ap*WO + WO*A + Np*WO*N + Cp*C = 0   ');
switch lower (reduce.balance.method)
    case 'iter'
        WO = math.glyap1(A,C,N,true,reduce.glyap);
    case 'bicg'
        WO = math.glyap2(A,C,N,true,reduce.glyap);
end

%% Determine method of balancing transform
switch lower (reduce.balance.transform)
    case 'srbt' % Square Root Balanced Truncation  
        [balanced.S, balanced.T, Sigma] = oct.srbt(WC, WO); 
    case 'mrmr' % Minima Realization and Model Reduction
        [balanced.S, balanced.T, Sigma] = oct.mrmr(WC, WO);
    otherwise
        util.error (['Wrong choice of balancing transform method : ' reduce.balance.transform])
end

% Plot details of balancing transformation
plot.balance (Sigma,WC,WO)

% Add extra blocks to S, T matrices to compensate for the splitting
% Apparently, this leaves T to be right-inverse of S, but not left-inverse 
if strcmpi(reduce.balance.A_stable,'ssu')
    dimt=size(balanced.T);
    balanced.S=[eye(n1) zeros(n1,dimt(1)); zeros(dimt(2),n1) balanced.S];
    balanced.S=balanced.S/U;
    balanced.T=[eye(n1) zeros(n1,dimt(2)); zeros(dimt(1),n1) balanced.T];
    balanced.T=U*balanced.T;
end

%% Carry out balancing transformation
% Check whether T is the (pseudo-)inverse of S 
util.pseudo(balanced.S,balanced.T)

% Transform A,B,N,C
balanced.A=balanced.S*bilinear.A*balanced.T;
for d=1:length(bilinear.B)
    balanced.B{d}=balanced.S*bilinear.B{d};
end
for d=1:length(bilinear.N)
    balanced.N{d}=balanced.S*bilinear.N{d}*balanced.T;
end
for d=1:length(bilinear.C)
    balanced.C{d}=bilinear.C{d}*balanced.T;
    balanced.Q{d}=bilinear.Q{d};
end

% State vectors (x, transformed) and output (y, unchanged)
balanced.x.initial=balanced.S*bilinear.x.initial;
balanced.x.equilib=balanced.S*bilinear.x.equilib;
balanced.y.initial=           bilinear.y.initial;
balanced.y.equilib=           bilinear.y.equilib;

% Labels of control targets and title
balanced.label=bilinear.label;  
balanced.title=bilinear.title;

% Save reduced matrices to data file
util.disp(['Saving balanced matrices A, B, N, C and vectors x,y to file : ' filename '_b.mat'])
save([filename '_b'],'balanced');

% Plot spectrum of A 
oct.spec_A (mfilename)

% Output clock/date/time
util.clock;

end
