% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '******************************************************************' )
util.disp ( 'Gaussian packet scattered from linear ramp' )
util.disp ( '******************************************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 192;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  20.0;              % Upper bound of grid

% Temporal discretization and propagator
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 020;                   % Index of final time step

time.main.delta  = 0.5;                  % Size of time steps 
time.sub.n       = 050;                  % Number of sub steps per time step

time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = -10.0;          % Lower truncation of energy
hamilt.truncate.max    = +30.0;          % Upper truncation of energy

hamilt.pot.handle      = @pot.taylor;    % Taylor series: Linear potential
hamilt.pot.v{1,1} = 1;                   % Slope parameter

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % Negative imaginary potential
hamilt.nip.exp  = 4;                     % Exponent
hamilt.nip.min = -07;                    % Beginning of inner-grid region
hamilt.nip.max =  15;                    % End of inner-grid region

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width = +1.0;                 % Width 
psi.dof{1}.pos_0 = -2.0;                 % Center in position representation
psi.dof{1}.mom_0 = +4.0;                 % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type    = 'contour';       % Contour plot of Wigner transform

plots.expect.energies.max = 20;          % manually set ranges for energy plot
