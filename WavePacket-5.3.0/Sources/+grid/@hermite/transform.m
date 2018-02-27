%------------------------------------------------------------------------------
%
% Transform the original, which is assumed to have the same number of points in
% the direction of this degree of freedom as the grid object obj, to a different
% number of points, keeping all the other variables constant.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function [ trafo, weights, xgrid ] = transform(obj, n_pts, original)

global space
persistent tmpdof tmpw tmpx

% The essential transformation is very straightforward.
%
% 1. We transform the original to FBR. The l-th component in the direction of
%    this DOF is then the coefficient for the l-th polynomial.
% 2. We extend the matrix to more or less polynomials. This can be done
%    by truncating it or by padding it with zeros in this direction.
% 3. We transform the result back to the DVR basis. Internally, this is done by
%    constructing a new temporary DOF object, and using it to transform the
%    result back. Since the fbr2dvr method does not use any global variables
%    except space.size.n_dim (which is the same for the transformed grid), this
%    is perfectly legal here.

% Throw out the trivial case
if n_pts == obj.n_pts
    trafo = original;
    weights = space.dvr.weight_1D{obj.dof};
    xgrid = space.dvr.grid_1D{obj.dof};
    return;
end

% Make a check that the input grid is ok
if ndims(original) ~= ndims(space.dvr.grid_ND{1})
    util.error('Could not transform grid because dimensions did not match')
end

%% 1. Transform to FBR

trafo = dvr2fbr(obj, original);


%% 2. Truncate or pad the grid in the direction of the current degree of
%%    freedom.

if n_pts > obj.n_pts
    % pad the array
    dimensions = size(trafo);
    dimensions(obj.dof) = n_pts - obj.n_pts;
    padding = zeros(dimensions);

    trafo = cat(obj.dof, trafo, padding);
else
    % truncate the array; We first bring it to the pseudo-2D form, truncate the
    % first dimension using Matlab functionality, and then perform some magic
    % to transform it back.
    [trafo, permutation, dimensions] = shape(obj, trafo);

    trafo = trafo(1:n_pts, :);

    dimensions(1) = n_pts;
    trafo = shape_back(obj, trafo, permutation, dimensions);
end


%% 3. Transform the result back to the DVR.

% Since the creation of a new grid object involves an expensive diagonalisation,
% and we typically do a lot of transformations in a row, we try to recycle this
% temporary object if possible.
if isempty(tmpdof) || tmpdof.n_pts ~= n_pts
    tmpdof = obj;
    tmpdof.n_pts = n_pts;
    [tmpdof, tmpw, tmpx, p] = init(tmpdof);
end

% Now transform back
trafo = fbr2dvr(tmpdof, trafo);
xgrid = tmpx;
weights = tmpw;
