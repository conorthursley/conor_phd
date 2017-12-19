% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.


function fbr = dvr2fbr(obj, dvr)

% Since matlab does not know how to handle matrix multiplication for more than
% 2 dimensions, we have to cast and recast the matrix to a pseudo-2D form.
% Apart from that, the transformation is just a matrix multiplication with a
% transformation matrix generated in grid_make.

[fbr, permutation, shapedims] = shape(obj, dvr);

fbr = obj.trafo2fbr * fbr;

fbr = shape_back(fbr, permutation, shapedims);
