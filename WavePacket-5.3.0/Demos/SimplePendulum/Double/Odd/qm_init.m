% Copyright (C) 2009 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '********************************************************' )
util.disp ( 'Evolution of double well pendulum with V = 100 -> 0     ' )
util.disp ( '********************************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % using fft grid
space.dof{1}.mass  = 1;                  % (Reduced) moment of inertia
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid 
space.dof{1}.x_max = 2*pi;               % Upper bound of grid

space.amo{1}.handle   = @amo.cosine;     % Projecting on cos-function
space.amo{1}.exp = 2;                    % Degree of alignment

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 200;                   % Index of final time step

time.main.delta = pi/100;                % Size of time steps
time.sub.n      =  05;                   % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min  = 000;              % Lower truncation of energy
hamilt.truncate.max  = 100;              % Upper truncation of energy

% Initial wave function
psi.dof{1}.handle   = @wav.pendulum;     % Mathieu type wavefunction
psi.dof{1}.parity   = 's';               % Sine elliptic (ground doublet)
psi.dof{1}.order    = 1;                 % Order of Mathieu function 
psi.dof{1}.multiple = 2;                 % Potential multiplicity
psi.dof{1}.barrier  = 100;               % Potential barrier
psi.dof{1}.shift    = 0;                 % Potential shift (horizontal)

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'polar';     % Contour plot of Wigner transform

plots.expect.energies.max = 4;           % Range for energy plot

plots.density.export.on = true;
plots.expect.export.on = true;
