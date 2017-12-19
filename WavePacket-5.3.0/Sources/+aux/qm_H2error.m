%--------------------------------------------------------------------------
%
% Calculated H2error from 'balanced truncation' and/or 'H2 model reduction'
% for equations of motions as determined by input variable "eom". The
% input vectors "truncate" and/or "reduce" contain a variable number of 
% approximation degrees to be simulated; if one of those vectors is left 
% empty on input, the respective simulation method will be skipped.
% Graphical output will be displayed in figure specified by integer "fig".
%
% switchoff=0: Do complete simulations
% switchoff>0: Skip calculation of bound states and matrix elements 
% switchoff>1: Skip calculation of A, B, N, C system matrices
% switchoff>2: Skip balancing transformation (BT method only)
%--------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2014    Tobias Breiten, TU Graz
% Copyright (C) 2015-16 Burkhard Schmidt, FU Berlin

function qm_H2error (eom, fig, truncate, singular, H2reduce, switchoff)

global reduce

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

% Set all parameters to default values if not specified by the user 
if ~isfield (reduce,'H2error')
    reduce.H2error=[];
end

if ~isfield(reduce.H2error,'A_shift')
    reduce.H2error.A_shift = 1e-6;
end
util.disp (['Value of A-shift : ' num2str(reduce.H2error.A_shift)])

if ~isfield(reduce.H2error,'BN_scale')
    reduce.H2error.BN_scale = 1;
end
util.disp(['Scaling B-vector, N-matrix: ' num2str(reduce.H2error.BN_scale)])

if ~isfield (reduce.H2error,'method')
    reduce.H2error.method='iter';
end
util.disp(['Method for GLYAP solving: ' reduce.H2error.method])


% Bound states, matrix elements (full dimensionality)
if switchoff<1
    qm_setup; qm_init; qm_bound;  qm_cleanup;
    qm_setup; qm_init; qm_matrix; qm_cleanup;
end

% ABNCD matrices (full dimensionality, ==> structure "bilinear")
if switchoff<2 
    qm_setup; qm_init; qm_abncd (eom); qm_cleanup;
else % if ABNCD already available
    load (eom)
end

%% Retrieve system matrices (full dimensionality, unbalanced)
global bilinear
A = sparse(bilinear.A);
n = size(A,1);

B = cell(length(bilinear.B),1);
for d=1:length(bilinear.B)
    B{d} = sparse(bilinear.B{d});
end

N = cell(length(bilinear.N),1);
for d=1:length(bilinear.N)
    N{d} = sparse(bilinear.N{d});
end

C = cell(length(bilinear.C),1);
for d=1:length(bilinear.C)
    C{d}=sparse(bilinear.C{d});
end
    
%% Simple truncation
if ~isempty (truncate)
    
    % Balancing transformation
    if switchoff<3
        qm_init; qm_balance (eom); qm_cleanup; 
        switchoff = 3;
    end
        
    % Loop: different truncations
    errBT = zeros(length(truncate),1);
    for k=1:length(truncate)
        qm_setup(); qm_init; qm_truncate(eom,'t',truncate(k));
        util.disp('   ')
        util.disp('--------------------------------')
        util.disp('Calculate error norm (BT method)')
        util.disp('--------------------------------')
        util.disp('   ')
        errBT(k) = err(A,B,N,C,n);
    end

    % Output clock/date/time
    util.clock;

end
%% Singular perturbation
if ~isempty (singular)
    
    % Balancing transformation
    if switchoff<3
        qm_init; qm_balance (eom); qm_cleanup; 
    end
        
    % Loop: different truncations
    errSP = zeros(length(singular),1);
    for k=1:length(singular)
        qm_setup(); qm_init; qm_truncate(eom,'s',singular(k));
        util.disp('   ')
        util.disp('--------------------------------')
        util.disp('Calculate error norm (SP method)')
        util.disp('--------------------------------')
        util.disp('   ')
        errSP(k) = err(A,B,N,C,n);
    end

    % Output clock/date/time
    util.clock;

end

%% H2 optimal model reduction
if ~isempty (H2reduce)
    
    % Loop: different reductions
    errH2 = zeros(length(H2reduce),1);
    for k=1:length(H2reduce)
        qm_setup(); qm_init; qm_H2model(eom,H2reduce(k));
        util.disp('   ')
        util.disp('--------------------------------')
        util.disp('Calculate error norm (H2 method)')
        util.disp('--------------------------------')
        util.disp('   ')
        errH2(k) = err(A,B,N,C,n);
    end
    
    % Output clock/date/time
    util.clock;
    
end

%% Initialize figure
global plots info reduce
figure (fig); clf
info.program = 'qm_H2error';
plot.logo;
mylegend = {};

% Plot BT result
if ~isempty (truncate)
    semilogy(truncate, real(errBT),'o','MarkerEdgeColor',plots.style.colors(1,:))
    mylegend{length(mylegend)+1} = 'Balanced truncation';
    hold on
end

% Plot BT result
if ~isempty (singular)
    semilogy(singular, real(errSP),'o','MarkerEdgeColor',plots.style.colors(2,:))
    mylegend{length(mylegend)+1} = 'Singular perturbation';
    hold on
end

% Plot H2 result
if ~isempty (reduce)
    semilogy(reduce, real(errH2),'o','MarkerEdgeColor',plots.style.colors(3,:))
    mylegend{length(mylegend)+1} = 'Interpolation based';
end

% General plot settings
legend(mylegend)
set ( gca, 'LineWidth',     plots.style.line.thick, ...
    'FontName',      plots.style.font.name,  ...
    'FontSize',      plots.style.font.large, ...
    'FontWeight',    plots.style.font.heavy )
xlabel (['dimensionality, reduced from ' int2str(n)])
ylabel ( 'H2 error' )
    
title  (['Error system || A shift: ' num2str(reduce.H2error.A_shift) ', B/N scaling: ' num2str(reduce.H2error.BN_scale) ', Glyap: ',reduce.H2error.method])

% Save figure and variables
saveas (gca, int2str(fig), 'fig')
save (int2str(fig))

end

function [H2err] = err (A,B,N,C,n)

% Retrieve system matrices (full dimensionality)
global bilinear reduce
Ar = sparse(bilinear.A);
r = size(Ar,1);

Br = cell(length(bilinear.B),1);
for d=1:length(bilinear.B)
    Br{d} = sparse(bilinear.B{d});
end

Nr = cell(length(bilinear.N),1);
for d=1:length(bilinear.N)
    Nr{d} = sparse(bilinear.N{d});
end

Cr = cell(length(bilinear.C),1);
for d=1:length(bilinear.C)
    Cr{d}=sparse(bilinear.C{d});
end

% Set up error system
AA = blkdiag(A, Ar);
AA = AA - reduce.H2error.A_shift * speye(n+r); % so far, only EVS implemented

BB = cell(length(B),1);
for d=1:length(B)
    BB{d} = [B{d}; Br{d}]/reduce.H2error.BN_scale;
end

NN = cell(length(N),1);
for d=1:length(N)
    NN{d} = blkdiag(N{d}, Nr{d})/reduce.H2error.BN_scale;
end

CC = cell(length(C),1);
for d=1:length(C)
    CC{d} = [C{d} -Cr{d}];
end

% Solve generalized Lyapunov equation for error system
switch reduce.H2error.method
    case 'iter'
        X = math.glyap1(AA,BB,NN,false,[]);
    case 'bicg'
        X = math.glyap2(AA,BB,NN,false,[]);
    otherwise
        util.error (['Wrong choice of GLYAP solver method : ' reduce.H2error.method])
end

% Calculate H2 error
H2err = 0;
for d=1:length(CC)
    H2err = H2err + trace(CC{d}*X*CC{d}');
end
H2err = sqrt(H2err);
util.disp(['H2 error norm = ', num2str(real(H2err))])
util.disp('   ')

end
