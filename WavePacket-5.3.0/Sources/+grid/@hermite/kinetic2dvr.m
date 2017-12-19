% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function dvrkin = kinetic2dvr(obj)

global space hamilt


if obj.nokin
    dvrkin = zeros(space.size.n_tot,space.size.n_tot);
    return;
end

if isempty(obj.kin)
    obj = init_kin(obj, 1);
end

% The situation here is very similar to the FFT one. We already have
% the kinetic energy given in DVR, and just need to upgrade it to a full
% multidimensional array. Essentially, the code there can be copied.

kinetic = obj.kin;

kinetic( abs(kinetic) < hamilt.eigen.cutoff ) = 0;
if hamilt.eigen.storage == 's'
    kinetic = sparse(kinetic);
end

% Transform to a Pseudo-2D matrix form
if hamilt.eigen.storage == 's'
    dvrkin = sparse(space.size.n_tot, space.size.n_tot);
else
    dvrkin = zeros(space.size.n_tot);
end

indices = cell(space.size.n_dim);
[indices{:}] = ind2sub(size(space.dvr.grid_ND{1}), (1:space.size.n_tot));

% Fill it row by row
for ii = 1:space.size.n_tot
    allowed    = ones(space.size.n_tot, 1);

    % Now disallow all indices, where the other dimension's indices change,
    % which translates to a zero matrix element
    for k = 1:space.size.n_dim
        if k == obj.dof
            continue;
        end

        allowed(indices{k} ~= indices{k}(ii)) = 0;
    end

    % get the indices of the allowed cells, and fill the row
    allowedindex = find(allowed);
    dvrkin(ii, allowedindex) = kinetic( indices{obj.dof}(ii), indices{obj.dof}(allowedindex) );
end
