%--------------------------------------------------------------------------
% Right-hand-side of bilinear control equation for state vector x(t):
%
%  d                     m                         
%  -- x(t) = A x(t) + i Sum u (t) ( N x(t) + b )           
%  dt                   k=1   k       k        k
%
% where A is used to describe the dynamics of the unperturbed system and 
% where the external control field(s) u_k(t) interact through N_k and b_k
%
% Here the time-dependence of u(t) is realized by calling util.efield(t)
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung, Burkhard Schmidt, Ulf Lorenz
%               2012 Burkhard Schmidt, Jeremy Rodriguez
%
% see the README file for license details.

function x_dot = rhs2 (t, x)

global bilinear

% dynamics of the unperturbed system
x_dot = bilinear.A*x;

% Get external control field
u = util.efield(t);

% interaction with external control fields along x and or y
if length(bilinear.N)>0 
    x_dot = x_dot + 1i * u.x  * ( bilinear.N{1} * x  + bilinear.B{1});
end
if length(bilinear.N)>1
    x_dot = x_dot + 1i * u.y  * ( bilinear.N{2} * x  + bilinear.B{2});
end

