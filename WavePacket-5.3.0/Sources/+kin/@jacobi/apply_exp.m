% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

function apply_exp(obj)

global psi space hamilt

% Basically the same as kin_apply, but with the short-time propagator.

for m = 1:hamilt.coupling.n_eqs
    % DVR => FBR
    fbr = dvr2fbr(space.dof{obj.dof_c}, psi.dvr.grid_ND{m});

    % Application of the grid and back-conversion
    psi.dvr.grid_ND{m} = fbr2dvr( space.dof{obj.dof_c}, fbr .* obj.grid_exp );
end
