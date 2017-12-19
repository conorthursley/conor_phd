% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = apply(obj, new)

% We do three things:
%
% 1. Convert the wavefunction to FBR
% 2. Multiply with our kinetic energy grid
% 3. Convert it back to DVR
%
% If new is set to true, we apply everything to psi.dvr.new_ND, otherwise
% to psi.dvr.grid_ND.
%
% Note that we transparently call kin_init if the grid is not yet set,
% so you might have to catch the initialised object in the output.

global psi space hamilt

% If kin_apply is called without previous initialisation, initialise it here.
if isempty(obj.grid)
	obj = init(obj, 1);
end


%% 1. Convert Wavefunction to FBR
for m = 1:hamilt.coupling.n_eqs
    if new
        psi.dvr.new_ND{m} = dvr2fbr(space.dof{obj.dof(1)}, psi.dvr.new_ND{m});
        psi.dvr.new_ND{m} = dvr2fbr(space.dof{obj.dof(2)}, psi.dvr.new_ND{m});
    else
        psi.dvr.new_ND{m} = dvr2fbr(space.dof{obj.dof(1)}, psi.dvr.grid_ND{m});
        psi.dvr.new_ND{m} = dvr2fbr(space.dof{obj.dof(2)}, psi.dvr.new_ND{m});
    end
end

%% 2. Apply the grid
for m = 1:hamilt.coupling.n_eqs
    psi.dvr.new_ND{m} = psi.dvr.new_ND{m} .* obj.grid;
end

%% 3. Convert back
for m = 1:hamilt.coupling.n_eqs
    psi.dvr.new_ND{m} = fbr2dvr(space.dof{obj.dof(1)}, psi.dvr.new_ND{m});
    psi.dvr.new_ND{m} = fbr2dvr(space.dof{obj.dof(2)}, psi.dvr.new_ND{m});
end
