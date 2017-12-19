% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Ulf Lorenz
%               2012 Ulf Lorenz
%
% see the README file for license details.

function dvrkin = kinetic2dvr(obj)

global space

% Note that the kinetic energy grid has already been assembled in k-space in
% the init() function, and is available under obj.grid. We only need to
% a) make a proper (diagonal) tensor out of it
% b) transform this tensor to the DVR using the matrix2dvr() method of the FFT
%    grids.
% c) produce a 2D matrix out of the higher-dimensional tensor.

% (a) make the grid a diagonal tensor with proper reshaping
req_size = size(obj.grid);
fbrkin = diag(obj.grid(:));
fbrkin = reshape(fbrkin, cat(2, req_size, req_size));

% b) transform it along both relevant degrees of freedom.
dvrkin = matrix2dvr(space.dof{obj.dof(1)}, fbrkin);
dvrkin = matrix2dvr(space.dof{obj.dof(2)}, dvrkin);

% c) the final output must be a 2D matrix, so reshape once again...
dvrkin = reshape(dvrkin, space.size.n_tot, space.size.n_tot);
