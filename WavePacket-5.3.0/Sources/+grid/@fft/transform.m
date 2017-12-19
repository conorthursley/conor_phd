%------------------------------------------------------------------------------
%
% Transform the original, which is assumed to have the same number of points in
% the direction of this degree of freedom as the grid object obj, to a different
% number of points, keeping all the other variables constant. Note that we do
% not change the boundaries of the grid!
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function [ trafo, weights, xgrid ] = transform(obj, n_pts, original)

global space

% The basic procedure is always the same.
% 
% 1. We first transform the original to the FBR basis
% 2. We add or truncate some coefficients
% 3. We transform it back.
%
% However, in contrast to, e.g., Legendre polynomials the truncation or
% expansion of the coefficients has to be done at the beginning and at the end
% of the array (i.e., for positive and negative k values).

% Throw out the trivial case
if n_pts == obj.n_pts
    trafo = original;
    weights = [];
    xgrid = [];
    return;
end

if round(n_pts/2) ~= n_pts/2
    util.error ('Do not want to handle odd grid points.')
end

%% 1. Transformation to FBR.
trafo = fftshift ( fft  ( original, [], obj.dof ), obj.dof );

%% 2. Addition or truncation of coefficients.

% Some general notes: The momenta are given by p_n = n * dp, where dp is the
% spectral resolution, and constant for both grids. n goes from -N/2 to N/2 - 2.
% So if we want to transform from grid points N1 to N2, we have to truncate/add
% |N2-N1|/2 points at the beginning and the end of the momentum grid.

if n_pts > obj.n_pts
    % pad with zeros
    gridsize = size(original);
    gridsize(obj.dof) = (n_pts - obj.n_pts) / 2;

    padding = zeros(gridsize);

    trafo = cat(obj.dof, padding, trafo, padding);
else
    % Stupid matlab that cannot treat ND arrays
    if space.size.n_dim > 1
        % truncate. Here we have the problem that I don't know of a generic matlab
        % mechanism to truncate a matrix along dimension k, so we recast it in a
        % pseudo-2D form (essentially what the shape function of grid.legendre does)
        permutorder = (1:space.size.n_dim);
        permutorder(obj.dof) = 1;
        permutorder(1) = obj.dof;
        trafo = permute(trafo, permutorder);

        dimensions = size(trafo);
        trafo = reshape(trafo, dimensions(1), prod(dimensions(2:end)));
    end

    % now truncate
    diff = (obj.n_pts - n_pts) / 2;
    trafo = trafo (diff+1:end-diff , :);

    if space.size.n_dim > 1
        % revert the casting and permutation
        dimensions(1) = n_pts;
        trafo = reshape(trafo, dimensions);
        trafo = ipermute(trafo, permutorder);
    end
end

%% 3. Scale the wave function. The scaling factor new_points/old_points
% needs to be applied because the IFFT divides by the (in this case incorrect)
% number of points.
trafo = trafo * (n_pts / obj.n_pts);
    
%% 4. Transform it back. We need to do FFT by hand because fbr2dvr
% is a bit picky about the grid sizes.

trafo = ifft ( ifftshift ( trafo, obj.dof ), [], obj.dof );
weights = ones(n_pts, 1) * (obj.x_max - obj.x_min) / n_pts;
xgrid = linspace(obj.x_min, obj.x_max - weights(1), n_pts)';
