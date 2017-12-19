% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function [shapedmat, permutation, shapedims] = shape(obj, input)

global space

% This function brings an input of the form of psi.dvr.grid_ND into a
% Pseudo-2D form needed because matlab cannot perform matrix
% multiplications for multi-dimensional arrays.

% We permute our original matrix such that the angular DOF becomes the
% first, and the radial DOF becomes the second matrix dimension.

if space.size.n_dim == 1
    % do nothing loop
    permutation = 1;
    shapedims = 1;
    shapedmat = input;
else
    % Get the permutation array
    permutation = (1:space.size.n_dim);
    permutation(1) = obj.dof;
    permutation(obj.dof) = 1;
    if ~isempty(obj.R_dof)
        permutation(2) = obj.R_dof;
        permutation(obj.R_dof) = 2;
    end

    % permute the matrix dimensions
    shapedmat = permute(input, permutation);

    % recast the matrix in pseudo-2D form
    shapedims = size(shapedmat);
    shapedmat = reshape(shapedmat, shapedims(1), prod(shapedims(2:end)));
end
