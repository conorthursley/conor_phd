% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%               based on the analogous function in grid.legendre
%
% see the README file for license details.

function [weights, xi] = quadrature (num_points)

% This function does the Gauss-Hermite quadrature and
% returns the weights, and quadrature points
% used for a DVR with an expansion in Harmonic oscillator eigenstates.
% Note that we use scaled coordinates, i.e., we assume the basis
% functions are exp(-y^2/2)*H_n(y). In real harmonic oscillator situations,
% the coordinate is y = sqrt(m \omega) * (x - x_0). This leads to a shift
% and a scaling of quadrature points and weights, though we put this
% into the grid_make function.

%% Find the grid points and weights. We employ the procedure proposed in
%% David Tannor's book.

ind = (1:num_points-1);

maindiag = zeros(num_points,1);
offdiag = sqrt(ind/2);

X = diag(maindiag) + diag(offdiag,-1) + diag(offdiag,1);
[eigvec,eigval] = eig(X);
xi = diag(eigval);
weights = (eigvec(1, : ).^2 / norm(eigvec(1, :)) * sqrt(pi))';
