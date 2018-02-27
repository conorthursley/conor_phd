%--------------------------------------------------------------------------
%
% Create initial state vector for TDSE propagation 
%
% Note that quantum states start counting from zero!
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%               2011 Boris Schaefer-Bung
%               2012 Jeremy Rodriguez, Ulf Lorenz
%
% see the README file for license details.
 
function tdse (H0)
global control bilinear

% Size of density matrices
dim=size(H0, 1);

util.disp (' ')
util.disp ('-------------------------------------------------------------')
util.disp (' ')

%% Initial state vector
psi = zeros(dim,1);
switch lower (control.initial.choice)
    case 'pure' % pure state
        ii = control.initial.pure+1;
        psi(ii) = 1;
        util.disp(['TDSE initial state: pure state = ', int2str(control.initial.pure)])
    case 'cat' % Schroedinger cat state
        ii = control.initial.cat(1)+1;
        jj = control.initial.cat(2)+1;
        psi(ii) = 1/sqrt(2);
        psi(jj) = 1/sqrt(2);
        util.disp(['TDSE initial state: cat state = ', int2str(control.initial.cat)])
    otherwise
        util.error (['Wrong choice of TDSE initial conditions: ' control.initial.choice])
end
bilinear.x.initial = psi;

%% Equilibrium state vector
psi    = zeros(dim,1);
psi(1) = 1;
bilinear.x.equilib = psi;
util.disp(['TDSE "equilibrium" : ground state'])
util.disp(' ')


