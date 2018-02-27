% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space

util.disp ( '*********************************************' )
util.disp ( 'Fulvene torsion: Stationary torsion (A state)' )
util.disp ( 'See http://dx.doi.org/10.1002/cphc.200600543 ' )
util.disp ( '*********************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1/(2.44/atomic.E.meV);  % Reduced moment of inertia: 1/2.44 meV
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -pi;                % Lower bound of grid 
space.dof{1}.x_max =  pi;                % Upper bound of grid

% Hamiltonian operator 
hamilt.eigen.symmetry  = 'g';            % Even parity eigenstates

hamilt.truncate.min    =  0.0;           % Lower truncation of energy
hamilt.truncate.max    = 0.25;           % Upper truncation of energy

hamilt.pot.handle    = @pot.pendulum;    % Intramolecular torsion
hamilt.pot.zeta = - 856*2.44/atomic.E.meV; % Barrier height: 856*2.44 meV

% Select eigen/values/functions
psi.eigen.start        = 0;              % Lower index
psi.eigen.stop         =10;              % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type          = 'contour'; % Make a contour plot of the eigenstates

plots.expect.energies.max = 0.05;        % manually set the ranges for the energy plot
