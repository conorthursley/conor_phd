% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2012 Ulf Lorenz
%
% see the README file for license details.

function obj = apply(obj, new)

% We do three things:
%
% 1. Convert the wavefunction to a mixed FBR representation for the
%    Legendre grid and a DVR one for the other grids.
% 2. Multiply with our kinetic energy grid. Since it is diagonal, this
%    is done pointwise instead of a real matrix multiplication.
% 3. Convert it back to pure DVR
%
% If new is set to true, we apply everything to psi.dvr.new_ND, otherwise
% to psi.dvr.grid_ND. The result is stored in psi.dvr.new_ND.
%
% Note that we transparently call kin_init if the grid is not yet set,
% so you might have to catch the initialised object in the output.

global psi space hamilt

% If kin_apply is called without previous initialisation, initialise it here.
if isempty(obj.grid)
    obj = init(obj, 1);
end

for m = 1:hamilt.coupling.n_eqs
    % DVR => FBR
    if new
        fbr = dvr2fbr(space.dof{obj.dof_c}, psi.dvr.new_ND{m});
    else
        fbr = dvr2fbr(space.dof{obj.dof_c}, psi.dvr.grid_ND{m});
    end

    % Apply the grid
    fbr = fbr .* obj.grid;

    % Convert it back.
    psi.dvr.new_ND{m} = fbr2dvr( space.dof{obj.dof_c}, fbr );
end
