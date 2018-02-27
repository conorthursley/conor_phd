% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2009 Ulf Lorenz
%               copied from grid.legendre
%
% see the README file for license details.

function dvr = matrix2dvr(obj, fbr)

% Copied from grid.legendre; should work here as well.

global space

% The matrix we want to transform has the shape
% A_{i1,i2,... j1,j2,...} = <i1,i2,...| A |j1,j2,...>
% where all indices denote the DVR basis except the one referring
% to this DOF's dimension, ik and jk, which is in FBR representation.
% Our task is to change this.

% This is done by applying the transformation two times, once for the
% ik index, once for jk. Each transformation consists of three steps:
% 1. Reshape the matrix A, to the form A_{ik/jk, r}, where r are all other
%    coordinates that we leave untouched.
% 2. Do the matrix multiplication with the appropriate transformation matrix
%    FBR=>DVR or its unitary transform
% 3. Reshape the matrix to the original form.

if length(size(fbr)) ~= 2*space.size.n_dim
    util.error('Input argument is not a matrix')
end

dvr = fbr;

%% Transform the first index.

% 1. Reshape the matrix. Similar to the things we do in shape.m
if space.size.n_dim > 1
    permutation = (1:2*space.size.n_dim);
    permutation(obj.dof) = 1;
    permutation(1) = obj.dof;

    dvr = permute(fbr, permutation);

    shapedims = size(dvr);
    dvr = reshape(dvr, shapedims(1), prod(shapedims(2:end)));
end

% 2. Multiply with the transformation matrix.
dvr = obj.trafo2dvr * dvr;

% 3. Shape it back to the original form
if space.size.n_dim > 1
    dvr = reshape(dvr, shapedims);
    dvr = permute(dvr, permutation);
end


%% Now transform the second index

% 1. Reshape again; Note that the second index transform the other way round,
%    so reshaping works a bit differently
if space.size.n_dim > 1
    permutation = (1:2*space.size.n_dim);
    permutation(obj.dof+space.size.n_dim) = 2*space.size.n_dim;
    permutation(end) = obj.dof + space.size.n_dim;

    dvr = permute(dvr, permutation);

    shapedims = size(dvr);
    dvr = reshape(dvr, prod(shapedims(1:end-1)), shapedims(end));
end

% 2. Multiply with the transformation matrix
dvr = dvr * obj.trafo2fbr;

% 3. Reshape back if necessary
if space.size.n_dim > 1
    dvr = reshape(dvr, shapedims);
    dvr = permute(dvr, permutation);
end
