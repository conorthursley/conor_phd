% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '********************************************************' )
util.disp ( 'Evolution of pendulum with squeezed initial function    ' )
util.disp ( '********************************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % using fft grid
space.dof{1}.mass  = 1;                  % (Reduced) moment of inertia
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid 
space.dof{1}.x_max = 2*pi;               % Upper bound of grid

space.amo{1}.handle   = @amo.cosine;     % Projecting on cos-function
space.amo{1}.exp = 1;                    % Degree of orientation

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 100;                   % Index of final time step

time.main.delta = pi/50;                 % Size of time steps
time.sub.n      =  05;                   % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min  = 000;              % Lower truncation of energy
hamilt.truncate.max  = 200;              % Upper truncation of energy

hamilt.pot.handle    = @pot.pendulum;    % Plane pendulum
hamilt.pot.eta = -50;                    % Prefactor of cos-potential
hamilt.pot.v_0 = +50;                    % Energy offset

% Initial wave function
psi.dof{1}.handle   = @wav.pendulum;     % Mathieu type wavefunction
psi.dof{1}.parity   = 'c';               % cosine elliptic  
psi.dof{1}.order    = 0;                 % Order of Mathieu function 
psi.dof{1}.multiple = 1;                 % Potential multiplicity
psi.dof{1}.barrier  = 100;               % Potential barrier
psi.dof{1}.shift    = pi/2;              % Potential shift (horizontal)

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'polar';     % Contour plot of Wigner transform
plots.expect.energies.max = 100;         % Range for energy plot
plots.density.export.on = true;
plots.expect.export.on = true;
