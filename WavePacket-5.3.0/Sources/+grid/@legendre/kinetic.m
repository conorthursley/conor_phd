% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function obj = kinetic(obj, usenew)

global psi hamilt

% The kinetic Hamiltonian is applied here in three steps
% 1. Transform the wave function to FBR
% 2. Multiply with a matrix representing the kinetic operator in FBR
% 3. Transform back the wave function to DVR
%
% Note that we use the internal shape and shape_back functions, this
% saves a few reshapes. The disadvantage is that obj.kin looks quite strange.

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
    % 1.
    if usenew == true
        [psi.dvr.new_ND{m}, permutation, shapedims] = shape(obj, psi.dvr.new_ND{m});
    else
        [psi.dvr.new_ND{m}, permutation, shapedims] = shape(obj, psi.dvr.grid_ND{m});
    end
    psi.dvr.new_ND{m} = obj.trafo2fbr * psi.dvr.new_ND{m};

    % 2.
    psi.dvr.new_ND{m} = obj.kin .* psi.dvr.new_ND{m};

    % 3.
    psi.dvr.new_ND{m} = obj.trafo2dvr * psi.dvr.new_ND{m};
    psi.dvr.new_ND{m} = shape_back(psi.dvr.new_ND{m}, permutation, shapedims);
end
