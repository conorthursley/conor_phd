% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function retval = momentum(obj, psi)

global space

if isempty(obj.kin)
    obj = init_kin(obj, 1);
end

% Applies the momentum operator to the input argument psi.
% Since this is usually done only for expectation values, i.e.
% not very often, we do not care about optimisation too much.

% Transform the wave function to FBR
retval = dvr2fbr(obj, psi);

% Apply the momentum. We assume it is stored in space.fbr.grid*
retval = retval .* space.fbr.grid_ND{obj.dof};

% Transform back
retval = fbr2dvr(obj, retval);
