% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.


function obj = kinetic(obj, usenew)

global psi hamilt

% The kinetic Hamiltonian is applied here in a single step, as the
% transformation to and from FBR has already been included in the
% definition of the internal matrix.

if isempty(obj.kin)
    obj = init_kin(obj, 1);
end

if obj.nokin
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.new_ND{m} = zeros(size(psi.dvr.grid_ND{m}));
    end
    return
end

for m = 1:hamilt.coupling.n_eqs
    % Shape the grid so that our DOF is first
    if usenew == true
        [psi.dvr.new_ND{m}, permutation, shapedims] = shape(obj, psi.dvr.new_ND{m});
    else
        [psi.dvr.new_ND{m}, permutation, shapedims] = shape(obj, psi.dvr.grid_ND{m});
    end

    % Apply the kinetic energy operator
    psi.dvr.new_ND{m} = obj.kin * psi.dvr.new_ND{m};

    % Reshape the grid to the original form
    psi.dvr.new_ND{m} = shape_back(psi.dvr.new_ND{m}, permutation, shapedims);
end
