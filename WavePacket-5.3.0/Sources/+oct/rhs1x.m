%--------------------------------------------------------------------------
% Right-hand-side of bilinear control equation for state vector x(t):
%
%  d                              (             )  
%  -- x(t) = A x(t) + i Sum u (t) ( N x(t) + b  )           
%  dt                    k   k    (  k        k )
%
% Note: x(t) is shifted with respect to equilibrium x_e with A x_e = 0 
% defining the equilibrium in the field-free case and with b = N x_e.
% Here A is used to describe the dynamics of the unperturbed system 
% and the external control field(s) u_k(t) interact through N_k and b_k
%
% Here the time-dependence of u(t) is realized by passing u(t) as argument
% 
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015 Burkhard Schmidt
%
% see the README file for license details.

function x_dot = rhs1x (u, x)

global bilinear

% dynamics of the unperturbed system
x_dot = bilinear.A*x;

% interaction with control fields u(t)
for d=1:length(u)
    x_dot = x_dot + 1i * u(d) *  ( bilinear.N{d} * x + bilinear.B{d} );
end

