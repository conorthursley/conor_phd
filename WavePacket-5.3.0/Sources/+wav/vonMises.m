%------------------------------------------------------------------------------
%
% This function creates the von-Mises distribution as a stationary state
% for combined orienting/aligning interaction of planar/spherical pendula
% see work published by B. Schmidt and B. Friedrich
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015 Burkhard Schmidt
%
% see the README file for license details.

function init_grid = vonMises (dir)
global space psi

% Some output
util.disp (' ')
util.disp ('*******************************************************')
util.disp ( ['Initial wavefunction for DOF :' int2str(dir)] )
util.disp ('   ' )
util.disp ('von Mises distribution')
util.disp ( ['Width parameter beta   : ' num2str(psi.dof{dir}.beta)] )

% Set up the grid
init_grid = exp(psi.dof{dir}.beta * cos(space.dvr.grid_ND{dir}));
