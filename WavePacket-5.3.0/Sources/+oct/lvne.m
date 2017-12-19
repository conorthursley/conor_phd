%--------------------------------------------------------------------------
%
% Create initial and thermal density matrix for LVNE propagation.
% Note that the matrices are vectorized: default is columnwise.
%
% Note that quantum states start counting from zero!
%
% For Boltzmann populations, temperature is in atomic units of 315,923.5 K
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
 
function lvne (H0)
global control bilinear

% Size of density matrices
dim=size(H0, 1);

util.disp (' ')
util.disp ('-------------------------------------------------------------')
util.disp (' ')

%% Initial density matrix
rho = zeros(dim);
switch lower (control.initial.choice)
    case 'pure' % pure state
        ii = control.initial.pure+1;
        rho(ii,ii) = 1;
        util.disp(['LvNE initial density: pure state = ', int2str(control.initial.pure)])
    case 'cat' % Schroedinger cat state: coherent superposition
        ii = control.initial.cat(1)+1;
        jj = control.initial.cat(2)+1;
        rho(ii,ii) = 1/2;
        rho(jj,jj) = 1/2;
        rho(ii,jj) = 1/2;
        rho(jj,ii) = 1/2;
        util.disp(['LvNE initial density: cat state = ', int2str(control.initial.cat)])
    case 'mixed' % Incoherent superposition
        ii = control.initial.mixed(1)+1;
        jj = control.initial.mixed(2)+1;
        rho(ii,ii) = 1/2;
        rho(jj,jj) = 1/2;
        util.disp(['LvNE initial density: mixed state = ', int2str(control.initial.mixed)])
    case 'thermal' % thermal (Boltzmann) distribution
        if (control.initial.temperature==0)
            rho(1,1) = 1;
        else
            boltz = exp(-H0/control.initial.temperature);
            rho = diag(boltz/sum(boltz));
        end
        util.disp(['LvNE initial density: thermal with kBT = ', num2str(control.initial.temperature)])
        for n=1:dim
            util.disp([int2str(n-1) ' : ' num2str(real(boltz(n)))])
        end
    otherwise
        util.error (['Wrong choice of LvNE initial conditions: ' control.initial.choice])
end
util.disp(' ')

% vectorize initial/thermal density matrices: columnwise ordering
bilinear.x.initial = rho(:);


%% Thermal density matrix as fixpoint of matrix A ==> Equilibrium
rho = zeros(dim);
if (control.lvne.temperature==0)
    rho(1,1) = 1;
else
    boltz = exp(-H0/control.lvne.temperature);
    rho = diag(boltz/sum(boltz));
end
util.disp(['LVNE equilibrium density: thermal with kBT = ', num2str(control.lvne.temperature)])
for n=1:dim
    util.disp([int2str(n-1) ' : ' num2str(rho(n,n))])
end
util.disp(' ')

% vectorize initial/thermal density matrices: columnwise ordering
bilinear.x.equilib = rho(:);

