% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it 
% into other works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '***********************************************' )
util.disp ( 'Razavy symmetric double well potential' )
util.disp ( 'Starting with Gaussian in the left well' )
util.disp ( 'Initial energy far below threshold' )
util.disp ( '***********************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 1/2;                 % Particle mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -7.0;               % Lower bound of grid 
space.dof{1}.x_max =  7.0;               % Upper bound of grid

% Hamiltonian operator 
hamilt.truncate.min    = -15.0;          % Lower truncation of energy
hamilt.truncate.max    = 150.0;          % Upper truncation of energy

hamilt.pot.handle      = @pot.razavy;    % Razavy hyperbolic potential
hamilt.pot.modified = true;
hamilt.pot.eta = -0.7;
hamilt.pot.zeta = 0.01;

% Temporal discretization and propagator
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 050;                   % Index of final time step

time.main.delta  = 0.2;                  % Size of time steps 
time.sub.n   =  100;                     % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Initial wave function (no initial momentum)
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width = 0.34;                 % Width of Gaussian
psi.dof{1}.pos_0 = -4.2483;              % Centered at left minimum (well)

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';
plots.density.pot.min = -15;
plots.density.pot.max = +15;
plots.density.contour.nlev = [30 30];
plots.expect.energies.max = 10;

