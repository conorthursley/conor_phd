% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function dvrkin = kinetic2dvr(obj)

global space hamilt


if obj.nokin
	dvrkin = zeros(space.size.n_tot,space.size.n_tot);
	return;
end


% use FGH method to analytically calculate the matrix elements
% See Tannor's book for details.
if obj.periodic == true
    T_nn = obj.p_max^2 / (6*obj.mass) * (1 + 2/obj.n_pts^2);
else
    T_nn = obj.p_max^2 / (6 * obj.mass);
end
kinetic = eye(obj.n_pts) * T_nn;

for n=1:obj.n_pts
    for m=1:n-1
        if obj.periodic == true
            T_nm = obj.p_max^2 * (-1)^(n-m) / (obj.mass * obj.n_pts^2 ...
                   * (sin((n-m)*pi/obj.n_pts))^2);
        else
            T_nm = obj.p_max^2*(-1)^(n-m) / ((pi*(n-m))^2 * obj.mass);
        end
        kinetic(n,m) = T_nm;
        kinetic(m,n) = T_nm;
    end
end

kinetic( abs(kinetic) < hamilt.eigen.cutoff ) = 0;
if hamilt.eigen.storage == 's'
    kinetic = sparse(kinetic);
end

% Now we basically have to upgrade our matrix T_ij to the generic form
% T_{i1,i2,...j1,j2,...} = eye(i1,i2) * eye(i2,j2) * ... * T_{ik,jk) * ...
% where the asterisk denotes the direct product. Unfortunately, MatLab
% is very poor with multidimensional arrays, and especially seems not to
% have anything like a direct product. Furthermore, if we want to allow
% sparse matrices, Matlab only handles them as two-dimensional objects, which
% prevents an "elegant" reshaping method to be used.
%
% As a basic solution, we setup the 2D matrix, then we go through each
% row, and fill it.

if hamilt.eigen.storage == 's'
    dvrkin = sparse(space.size.n_tot, space.size.n_tot);
else
    dvrkin = zeros(space.size.n_tot, space.size.n_tot);
end

indices = cell(space.size.n_dim, space.size.n_dim);
[indices{:}] = ind2sub(size(space.dvr.grid_ND{1}), (1:space.size.n_tot));

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
