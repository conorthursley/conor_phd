% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '*****************************' )
util.disp ( 'Morse oscillator (OH radical)' )
util.disp ( '*****************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1728.539;            % Reduced mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = 1.0;                % Lower bound of grid 
space.dof{1}.x_max = 10.0;               % Upper bound of grid
space.dof{1}.periodic = false;           % Build the kinetic energy matrix
                                         % without periodic boundary conditions

% Hamiltonian operator 
hamilt.truncate.min    =  0.0;           % Lower truncation of energy
hamilt.truncate.max    =  1.0;           % Upper truncation of energy

hamilt.pot.handle      = @pot.morse;     % Harmonic oscillator
hamilt.pot.d_e  = 0.1994;         % Dissociation energy
hamilt.pot.r_e  = 1.821;          % Equilibrium length
hamilt.pot.alf  = 1.189;          % Range parameter
hamilt.pot.t_e  = 0.0;            % Energetic shift

% Select eigen/values/functions
psi.eigen.start        = 00;             % Lower index
psi.eigen.stop         = 21;             % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'contour';   % Contour plot of Wigner function

plots.expect.energies.max = 0.25;        % Manually set range for energy plot
