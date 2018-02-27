%--------------------------------------------------------------------------
%
% Extract psi from columns(!) of eigenvector matrix and normalize
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function eigen ( step )
global hamilt psi space

index = step + psi.eigen.start;

% With symmetry
if isfield(hamilt.eigen, 'symmetry')
    bloated = hamilt.eigen.transform' * hamilt.eigen.eig_vecs(:, index);
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = bloated((m-1)*space.size.n_tot+1:m*space.size.n_tot);
    end

% General case: without symmetry
else
    for m = 1:hamilt.coupling.n_eqs
        psi.dvr.grid_ND{m} = hamilt.eigen.eig_vecs((m-1)*space.size.n_tot+1:m*space.size.n_tot,index);
    end
end

% Reshape eigenvector for >1 dimension
for k = 1:space.size.n_dim
    dims(k) = space.dof{k}.n_pts;
end
if space.size.n_dim > 1
    psi.dvr.grid_ND{m} = reshape(psi.dvr.grid_ND{m}, dims);
end


% Normalize eigenfunction
norm2 = 0;
for m = 1:hamilt.coupling.n_eqs
    norm2 = norm2 + sum(abs(psi.dvr.grid_ND{m}(:)).^2 .* space.dvr.weight_ND(:));
end

for m = 1:hamilt.coupling.n_eqs
    psi.dvr.grid_ND{m} = psi.dvr.grid_ND{m} / sqrt(norm2);
end
