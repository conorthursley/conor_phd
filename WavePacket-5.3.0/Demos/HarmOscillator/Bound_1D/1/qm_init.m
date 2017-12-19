% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '************************************' )
util.disp ( 'Harmonic oscillator in one dimension' )
util.disp ( '************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % Fourier grid
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = -7;                 % Lower bound of grid
space.dof{1}.x_max =  7;                 % Upper bound of grid
space.dof{1}.mass  =  1;                 % Mass for the kinetic energy

% Hamiltonian operator 
% hamilt.eigen.symmetry='u';

hamilt.truncate.min    =  0.0;           % Lower truncation of energy
hamilt.truncate.max    = 25.0;           % Upper truncation of energy

hamilt.pot.handle = @pot.taylor;         % Taylor expansion
hamilt.pot.v{1,1} = [0;1];        % Force constant

% Select eigen/values/functions
psi.eigen.start        =  0;             % Lower index
psi.eigen.stop         = 10;             % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.surface.view  = [30 60];   % View point for surface plot: [az el]

plots.density.type = 'contour';          % contour plot

plots.density.pot.max = 10;              % Dirty tweaking of density display
plots.expect.energies.max = 15;          % Tune the display of the energy values
