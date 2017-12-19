%------------------------------------------------------------------------------
%
% Applies the momentum operator to the input argument psi.
% Since this is done relatively rarely for calculating expectation values, 
% it is not too optimised.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007 Martin Winter
%               2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function retval = momentum(obj, psi)

global space

% Transform the wave function to FBR
retval = dvr2fbr(obj, psi);

% Apply the momentum operator: Note that space.fbr.grid_ND stores
% sqrt(l*(l+1))
retval = retval .* space.fbr.grid_ND{obj.dof};

% Transform back
retval = fbr2dvr(obj, retval);
