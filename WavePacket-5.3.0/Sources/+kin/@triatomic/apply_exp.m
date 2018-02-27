% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function apply_exp(obj)

% We do three things:
%
% 1. Convert the wavefunction to FBR
% 2. Multiply with our exponentiated kinetic energy grid
% 3. Convert it back to DVR
%
% If new is set to true, we apply everything to psi.dvr.new_ND, otherwise
% to psi.dvr.grid_ND.

global psi space hamilt


%% 1. Convert Wavefunction to FBR
for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = dvr2fbr(space.dof{obj.dof(1)}, psi.dvr.grid_ND{m});
        psi.dvr.grid_ND{m} = dvr2fbr(space.dof{obj.dof(2)}, psi.dvr.grid_ND{m});
end

%% 2. Apply the grid
for m = 1:hamilt.coupling.n_eqs
    psi.dvr.grid_ND{m} = psi.dvr.grid_ND{m} .* obj.grid_exp;
end

%% 3. Convert back
for m = 1:hamilt.coupling.n_eqs
    psi.dvr.grid_ND{m} = fbr2dvr(space.dof{obj.dof(1)}, psi.dvr.grid_ND{m});
    psi.dvr.grid_ND{m} = fbr2dvr(space.dof{obj.dof(2)}, psi.dvr.grid_ND{m});
end
