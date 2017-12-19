% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '*************************************' )
util.disp ( 'Harmonic oscillator in two dimensions' )
util.disp ( '*************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % Fourier grid
space.dof{1}.n_pts = 48;                 % Number of grid points
space.dof{1}.x_min = -7;                 % Lower bound of the grid
space.dof{1}.x_max =  7;                 % Upper bound of the grid
space.dof{1}.mass  =  1;                 % Mass for the kinetic energy

space.dof{2}       = grid.fft;           % Fourier grid
space.dof{2}.n_pts = 48;                 % Number of grid points
space.dof{2}.x_min = -7;                 % Lower bound of the grid
space.dof{2}.x_max =  7;                 % Upper bound of the grid
space.dof{2}.mass  =  1;                 % Mass for the kinetic energy

% Hamiltonian operator 
hamilt.eigen.storage   = 's';            % sparse storage

hamilt.truncate.min    =  0.0;           % Lower truncation of energy
hamilt.truncate.max    = 50.0;           % Upper truncation of energy


hamilt.pot.handle = @pot.taylor;         % Taylor expansion
hamilt.pot.v{1,1} = [0 0; 0.81 1];% Force constants

% Select eigen/values/functions
psi.eigen.start        = 000;            % Lower index
psi.eigen.stop         = 025;            % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.surface.view  = [25 55];   % View point for surface plot: [az el]
plots.density.type = 'surface';          % 3D surface plot
plots.expect.energies.max = 15;          % manually set range for energy plot
