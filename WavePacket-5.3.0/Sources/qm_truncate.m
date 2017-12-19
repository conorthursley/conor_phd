%------------------------------------------------------------------------------
%
% Truncate/reduce A, B, N, C and initial/equililibrium state vectors
%
% Upon balancing transformation, the obtained Hankel singular values (HSVs)
% indicate controllability and observability, and the resulting modes are 
% ordered accordingly. To reduce dimensionality, there are two methods: 
%
% Method = 't' (truncatation): 
% In simple tuncation one simply eliminates the low HSV modes because
% they are weakly controllable and weakly observable at the same time. 
%
% Method = 's' (singular perturbation):
% The idea of confinement truncation is based on the analogy 
% between large HSV-modes with slow dof's and low HSV-modes
% with fast dof's. Then an averaging principle ( based on singular 
% perturbation theory) can be used to derive equations of motion 
% for the former one, where the latter one is confined to its average, 
% or rather, the t --> oo limit.  
%
% A = A11 - A12 A22 \ A21
% N = A11 - N12 A22 \ A21
% C = C1  - C2  A22 \ A21
%
% From both an execution time and numerical accuracy standpoint, it is better
% to use Matlab's matrix division operator mldivide (A22\...) than inv(A22)
%
% Already balanced A, B, N, C matrices and initial/equilibrium state vectors 
% are read from "filename_b.mat" where filename typically is 'lvne' or 'tdse'. 
% In the end, the truncated matrices and vectors are written to either 
% "filename_t<n>.mat" (method='t') for (simple) truncation or 
% "filename_s<n>.mat" (method='s') for singular perturbation theory.
% Input parameter 'dim' specifies the dimension of the truncated system.
%
% C. Hartmann, B. Schaefer-Bung, A. Thöns-Zueva
% SIAM J. Control Optim., 51(3), 2356-2378 (2013)
% doi:10.1137/100796844
% 
% C. Hartmann, V. Vulcanov, Ch. Schütte
% Multiscale Model. Simul., 8(4), 1348-1367 (2010)
% doi:10.1137/080732717
% 
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung
% Copyright (C) 2012-16 Burkhard Schmidt
%
% see the README file for license details.
% 

function qm_truncate (filename, method, dim)

global  balanced bilinear

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

% Load full matrices (already balanced) from data file
load ([filename '_b']);

util.disp (' ')
util.disp ('-------------------------------------------------------------')
util.disp (' Truncating A, B, N, C matrices; also x_i, x_e state vectors ')
util.disp (' http://sourceforge.net/p/wavepacket/matlab/wiki/Reference.Programs.qm_truncate');
util.disp ('-------------------------------------------------------------')
util.disp (' ')
util.disp(['truncating to ' int2str(dim) ' modes:'])

nbal=size(balanced.A,1);

switch lower(method)
    case 's' % singular perturbation for A, N, C
        
        util.disp('singular perturbation approximation')
        bilinear.A=balanced.A(1:dim,1:dim)...            % A11
            - balanced.A(1:dim,dim+1:nbal)* ...          % A12
            ((balanced.A(dim+1:nbal,dim+1:nbal))...      % A22^(-1)
            \balanced.A(dim+1:nbal,1:dim));              % A21
        for d=1:length(balanced.N)
            bilinear.N{d}=balanced.N{d}(1:dim,1:dim)...  % N11
                -balanced.N{d}(1:dim,dim+1:nbal)*...     % N12
                ((balanced.A(dim+1:nbal,dim+1:nbal))...  % A22^(-1)
                \balanced.A(dim+1:nbal,1:dim));          % A21
        end
        for d=1:length(balanced.C)
            bilinear.C{d}=balanced.C{d}(1:dim)...           % C1
                -balanced.C{d}(dim+1:nbal)* ...             % C2
                ((balanced.A(dim+1:nbal,dim+1:nbal))...  % A22^(-1)
                \balanced.A(dim+1:nbal,1:dim));          % A21
            bilinear.Q{d}=balanced.Q{d};

        end
        
        
    case 't' % truncation  for A, N, C
        
        util.disp('simple truncation')
        bilinear.A=balanced.A(1:dim,1:dim);
        for d=1:length(balanced.N)
            bilinear.N{d}=balanced.N{d}(1:dim,1:dim);
        end
        for d=1:length(balanced.C)
            bilinear.C{d}=balanced.C{d}(1:dim);
            bilinear.Q{d}=balanced.Q{d};
        end
        
    otherwise
        util_error (['Invalid choice of truncation method : ' method])
end

oct.check_stable ('truncated system matrix A', bilinear.A)

% truncate/reduce B matrices
for d=1:length(balanced.B)
    bilinear.B{d}=balanced.B{d}(1:dim,:);
end

% truncate transformation matrices and state vectors (initial and equilibrium) 
bilinear.S=balanced.S(1:dim,:); % needed for qm_correlate ==> oct.reconstruct
bilinear.T=balanced.T(:,1:dim); % needed for qm_correlate ==> oct.reconstruct
bilinear.x.initial=balanced.x.initial(1:dim,:);
bilinear.x.equilib=balanced.x.equilib(1:dim,:);
bilinear.y.initial=balanced.y.initial;
bilinear.y.equilib=balanced.y.equilib;

% Save truncated/reduced matrices to data file
bilinear.label=balanced.label; % labels of control targets
switch lower(method)
    case 't'
        filename = [filename '_t' int2str(dim)];
        bilinear.title=[balanced.title ', simple truncation to '];
    case 's'
        filename = [filename '_s' int2str(dim)];
        bilinear.title=[balanced.title ', sing. pert. theory to '];
end

util.disp(['Saving truncated matrices A, B, N, C and densities to file : ' filename '.mat'])
save (filename, 'bilinear')

% Plot spectrum of A 
oct.spec_A (mfilename)

% Output clock/date/time
util.clock;

end

