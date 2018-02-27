% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function reshaped_mat = shape_back(input, permutation, shapedims)

global space

% Casts a matrix from the pseudo-2D form created by shape into the proper
% multidimensional form again. permutation and shapedims are passed from the
% original reshaping, so we do not have to duplicate the code here.

if space.size.n_dim == 1
    % Nothing to be done.
    reshaped_mat = input;
else
    % Shape the matrix back to ND form
    reshaped_mat = reshape(input, shapedims);

    % permute the matrix dimensions back
    reshaped_mat = ipermute(reshaped_mat, permutation);
end
