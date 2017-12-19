% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%
% see the README file for license details.

function kinetic_exp(obj)

% Performs a propagation step in the split operator method.
% Basically the same as grid_kinetic, but uses
% obj.kin_expo for the multiplication.

global psi hamilt

if obj.nokin
    return
end

% Fourier-transform, apply the factor, and transform back.
for m = 1:hamilt.coupling.n_eqs
    psi.dvr.grid_ND{m} = ifft( obj.kin_expo .* fft(psi.dvr.grid_ND{m}, [], obj.dof), [], obj.dof);
end
