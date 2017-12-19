% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function [weights, xi, trafo] = quadrature (num_points, m)

% This function does the Gauss-Jacobi quadrature and
% returns the weights, quadrature points, and transformation matrix
% used for a DVR with Spherical harmonics. The arguments are the number of
% points for the quadrature (NOT THE MAXIMUM ANGULAR MOMENTUM!), and the fixed
% value for the second quantum number.

%% Find the grid points and weights.

ind = (1:num_points-1);

% We use Gauss-Jacobi quadrature, which is the more general form for the
% associated Legendre "polynomials" with m ~= 0. Eigenvalues are the
% quadrature points, the weights are obtained from the first element of each
% eigenvector.  Note that we do not need the explicit form of the Jacobi
% polynomials, just their recursion relation.

maindiag = zeros(num_points,1);
offdiag = sqrt(ind.*(ind+2*m)./(2*ind+2*m-1)./(2*ind+2*m+1));

X = diag(maindiag) + diag(offdiag,-1) + diag(offdiag,1);
[eigvec,eigval] = eig(X);
xi = diag(eigval);
weights = eigvec(1, : ).^2' / norm(eigvec(1, :)) * sum(factorial(m)*...
            (-ones(1,m+1)).^(0:m)*2./(2*(0:m)+1) ./...
            factorial(0:m)./factorial(m-(0:m))) ./ ((1-xi.*xi).^(m));

%% Caclulate the transformation matrices.

% It can easily be shown that the transformation Psi(x_i) => f_l
% for an expansion in spherical harmonics Psi(x) = sum_l f_l Y_lm(x)
% is done by a matrix multiplication f_l = sum_i T_li * Psi(x_i) with
% T_li = Y_lm(x_i).

trafo = zeros(num_points, num_points);
for i = 1:num_points
    % Again, matlab confronts us with problems: legendre() returns all
    % values for _fixed_ l and all possible m, so we have to extract the value
    % of interest. Despite the name, the third argument returns spherical
    % harmonics up to a factor of sqrt(2pi).
    Ym_lm = legendre(i+m-1, xi, 'norm');
    trafo(i, :) = Ym_lm(m+1, :);
end

% Up to now, we only got the weights etc. for an expansion in normalized
% associated Legendre functions (/special Jacobi polynomials). However,
% we would like to operate with spherical harmonics. Since they have the
% second coordinate, phi, over which one has to integrate, the harmonics
% themselves are smaller by a factor of 1/sqrt(2pi), while the weight has to
% be increased by 2pi to reflect the additional integration over d_phi.
% This seems to be a purely academic shifting of factors, but the advantage is
% that the DVR then really tells us the wave function at the position, not the
% wave function scaled by sqrt(2pi), which can be slightly annoying in various
% circumstances.
weights = 2*pi * weights;
trafo   = trafo / sqrt(2*pi);
