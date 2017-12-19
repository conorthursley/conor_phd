% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function [ obj, weights, x_grid, p_grid ] = init_grid(obj)

% The stuff itself is calculated by Gaussian quadrature in an internal function.
[weights, x_grid, obj.trafo2fbr] = quadrature(obj.l_max - abs(obj.m_0) + 1, abs(obj.m_0));

% The spectral grid is sqrt( l(l+1) )
p_grid = (abs(obj.m_0):obj.l_max)';
p_grid = sqrt(p_grid.*(p_grid+1));

% The inverse trafo is just the transpose of the simple trafo, since it is
% unitary (even orthogonal)
obj.trafo2dvr = obj.trafo2fbr';

% Ah, yes. Don't forget to add the weights to the transformation matrix, since it
% involves an integration/summation over the DVR grid.
obj.trafo2fbr = obj.trafo2fbr .* repmat(weights', [ length(weights) 1] );
