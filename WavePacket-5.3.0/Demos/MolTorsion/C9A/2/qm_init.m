% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init(symmetry)
global atomic hamilt plots psi space

util.disp ( '*********************************************' )
util.disp ( 'Calculation of the absorption spectrum of the' )
util.disp ( 'S1 state of 9-(N-carbazolyl)-anthracene (C9A)' )
util.disp ( 'see Z.Phys.D 34:111                          ' )
util.disp ( '*********************************************' )

% moment of intertia is hbar/(4pi*c * 0.035 cm^-1)
c = 299792458/atomic.l.m*atomic.t.s; % approx 137
inertia = 1/(4 * pi * c * 0.035 * atomic.l.A * 1e-8);

% Spatial discretization
space.dof{1}       = grid.fft;           % Fourier grid
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = 030*pi/180;         % Lower bound of grid 30 degrees
space.dof{1}.x_max = 150*pi/180;         % Upper bound of grid 150 degrees
space.dof{1}.mass  = inertia;            % Moment of inertia for the kinetic energy

% Hamiltonian operator 
hamilt.eigen.symmetry  = symmetry;       % Take the symmetry from a global variable.

hamilt.truncate.min    = -1.0;           % Lower truncation of energy
hamilt.truncate.max    = 1.0;            % Upper truncation of energy

hamilt.pot.handle      = @pot.C9A;
hamilt.pot.state= 2;              % S1 state of C9A

% Select eigen/values/functions
psi.eigen.start        =  0;             % Lower index
psi.eigen.stop         =  35;            % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % contour plot

plots.density.on = false;                % Turn off uninteresting plotting
plots.expect.on  = false;
