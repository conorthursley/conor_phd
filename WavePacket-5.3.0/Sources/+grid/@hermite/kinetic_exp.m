% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function kinetic_exp(obj)

global psi hamilt

% Does the same as grid_kinetic, but uses an exponentiated form for the
% multiplication for the split operator method.

if obj.nokin
    return
end

for m = 1:hamilt.coupling.n_eqs
    % 1. Shape the grid to a suitable form
    [psi.dvr.grid_ND{m}, permutation, shapedims] = shape(obj, psi.dvr.grid_ND{m});

    % 2. Apply the prepared propagator
    psi.dvr.grid_ND{m} = obj.kin_expo * psi.dvr.grid_ND{m};

    % 3. Reshape the grid to the original form
    psi.dvr.grid_ND{m} = shape_back(obj, psi.dvr.grid_ND{m}, permutation, shapedims);
end
