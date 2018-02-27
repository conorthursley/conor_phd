% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%
% see the README file for license details.

function [ obj, weights, x_grid, p_grid ] = init_grid(obj)

%% The stuff itself is calculated by Gaussian quadrature in an internal function.
[weights, x_grid] = quadrature(obj.n_pts);


%% The frequency can be given either as force constant or as angular frequency.
%% synchronize the two. v_2 = mass * omega^2
if isempty(obj.omega)
    obj.omega = sqrt(obj.v_2 / obj.mass);
else
    obj.v_2 = obj.mass * obj.omega^2;
end

%% Toy around with the weights.

% We would like the weights in such a way that e.g., the norm can be
% calculated via sum_i w_i |\Psi(x_i)|^2.
% Unfortunately, this is not possible up to here. The reason is that
% the real expression N =  \int |\Psi(x)|^2 dx   cannot be converted to
% the sum directly , since \Psi is not a polynomial. So again we insert a
% one to get for the integral:
% N = \int exp(-x^2) | \Psi(x) exp(x^2/2) |^2 dx
% The expression in the absolute value is now a polynomial and can be expanded
% over the sum to give
% N = \sum_i w_i exp(x_i^2) |\Psi(x_i)|^2.
% As a result, to use the weights in a _convenient_ way, we have to include
% the additional factor exp(x_i^2).

weights = weights .* exp(x_grid.^2);


%% Scale coordinates

% Disadvantage: The Gaussian quadrature treats the special problem m*omega = 1,
% r_e = 0. To map this on unscaled coordinates, we have to do three things:

factor = sqrt(obj.mass * obj.omega);

% 1. Scale the coordinates x_real = x_orig / sqrt(m * omega)
x_grid = x_grid / factor;

% 2. Shift the coordinates to the real center of the potential
x_grid = x_grid + obj.r_e;

% 3. Scale the weights.
weights = weights / factor;

% Now we are (almost) ready.


%% For the spectral grid, we use the quantum numbers.
p_grid = (0:obj.n_pts-1)';


%% Construct transformation matrices

% The transformation matrices aim at providing means of converting
% DVR to FBR and vice versa (f is the scaling factor, N normalisation).
% DVR: \psi.i  = \Psi(x_i)
% FBR: \Psi(x) = \sum_n c_n 1/N exp(-(fx)^2/2) H_n(fx)
%
% DVR => FBR c_n = \int 1/N exp(-(fx)^2/2) H_n(fx) * \Psi(x) dx = \sum_i T_ni \psi.i
% FBR => DVR \psi.i = \sum_n c_n 1/N exp(-(fx)^2/2) H_n(fx) = \sum_n R_in c_n
%
% The back-conversion FBR=>DVR is very simple, the matrix is obviously
% R_in = 1/N exp(-(fx)^2/2) H_n(fx).
%
% The conversion FBR=>DVR can be converted to a sum over the quadrature
% points; with all the previous transformations, this works now.
% c_n = \sum_i w_i 1/N exp(-(fx)^2/2) H_n(fx) \psi.i
% so that
% T_ni = w_i 1/N exp(-(fx)^2/2) H_n(fx)
% Note that math.hermite already gives normalized Hermite polynomials for
% scaled coordinates (f = 1), so we only have to apply a factor sqrt(f).

obj.trafo2fbr = zeros(obj.n_pts);
obj.trafo2dvr = zeros(obj.n_pts);
for n = 1:obj.n_pts
    obj.trafo2fbr(n, :) = sqrt(factor) * weights .* exp(-(factor*(x_grid-obj.r_e)).^2/2)...
                             .* math.hermite(factor * (x_grid-obj.r_e), n-1);
    obj.trafo2dvr(:, n) = sqrt(factor) * exp(-(factor*(x_grid-obj.r_e)).^2/2) ...
                                  .* math.hermite(factor * (x_grid-obj.r_e), n-1);
end
