%------------------------------------------------------------------------------
%
% This function provides a very rudimentary and ugly functionality for
% transforming wave functions between different grid sizes.
%
% The background here is that one sometimes faces the problem that e.g.
% two calculations done with different grids have to be compared. Or that
% a calculation can be effieciently done with very few grid points, but
% the evaluation of an expectation value requires many more points. In these
% cases, we have to transform the wavefunction.
%
% The function takes three input arguments. "original" is the wave function
% to be transformed, which is assumed to fit to the grids under space.dvr.*.
% "newsize" is an array that gives the new number of points in each direction.
% The function then returns the transformed wavefunction ("trafo"). If 
% "outputgrid" is set to true, it also calculated the new grids and weights, and
% returns them as a second return value. Otherwise, this second value should
% be ignored.
%
% Limitations:
% - it is quite ugly
% - it cannot transform between different grids. I.e., if a degree of freedom
%   is an FFT grid, one cannot transform this to a Gauss-Hermite basis.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function [ trafo, newgrid ] = transform(original, newsize, outputgrid)

global space

%% Catch errors first
if size(newsize) ~= space.size.n_dim
    util.error ('Dimensions do not match. Cannot do transformation.')
end


%% Now transform the wave function.
trafo = original;
for k = 1:space.size.n_dim
    if newsize(k) == size(space.dvr.grid_1D{k}(:))
        newgrid.weight_1D{k} = space.dvr.weight_1D{k};
        newgrid.grid_1D{k} = space.dvr.grid_1D{k};
    else
        [trafo, newgrid.weight_1D{k}, newgrid.grid_1D{k}] = ...
                    transform(space.dof{k}, newsize(k), trafo);
    end
end

%% If output is requested, do some additional calculation for the newgrid stuff
if outputgrid
    % Calculate the multidimensional grids.
    if space.size.n_dim == 1
        newgrid.grid_ND{1} = newgrid.grid_1D{1};
        tmpweights{1} = newgrid.weight_1D{1};
    else
        [newgrid.grid_ND{1:space.size.n_dim}] = ndgrid(newgrid.grid_1D{1:space.size.n_dim});
        [tmpweights{1:space.size.n_dim}] = ndgrid(newgrid.weight_1D{1:space.size.n_dim});
    end

    % Calculate the total weight function.
    if space.size.n_dim == 1
        newgrid.weight_ND = tmpweights{1};
    else
        newgrid.weight_ND = tmpweights{1};
        for k = 2:space.size.n_dim
                newgrid.weight_ND = newgrid.weight_ND .* tmpweights{k};
        end
    end
end
