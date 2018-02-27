% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Ulf Lorenz
%               2012 Ulf Lorenz
%
% see the README file for license details.

function dvrkin = kinetic2dvr(obj)

global space

% Note that we already have a (grid) representation of the kinetic energy in
% the obj.grid variable (see apply.m). We only need to blow this up to a
% higher-order tensor, and transform it to the DVR along the radial degree
% of freedom.
%
% Note: the kinetic energy must have been initialized before calling this function.

%% 1. Make a matrix out of the kinetic energy from the internal grid
req_size = size(obj.grid);
fbrkin = diag(obj.grid(:));
fbrkin = reshape(fbrkin, cat(2, req_size, req_size));

%% 2. Transform it to the DVR
dvrkin = matrix2dvr(space.dof{obj.dof_c}, fbrkin);

%% 3. Reshape to make it a matrix
dvrkin = reshape(dvrkin, space.size.n_tot, space.size.n_tot);
