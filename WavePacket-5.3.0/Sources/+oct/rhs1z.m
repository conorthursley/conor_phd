%--------------------------------------------------------------------------
% Right-hand-side of bilinear control equation for Lagr. multiplier z(t):
%
%  d            +                     +               
%  -- z(t) = - A z(t) + i Sum u (t)  N z(t)        
%  dt                      k   k      k 
%
% where A is used to describe the dynamics of the unperturbed system 
% and the external control field(s) u_k(t) interact through N_k 
% Note: z(t) is *not* shifted with respect to equilibrium z_e with A^+z_e=0
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

function z_dot = rhs1z (u, z)

global bilinear

% dynamics of the unperturbed system
z_dot = - bilinear.A'*z;

% interaction with control fields u(t)
for d=1:length(u)
    z_dot = z_dot + 1i * u(d) * bilinear.N{d}' * z;
end

