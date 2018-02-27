% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '***************************************' )
util.disp ( 'Squeezed state of a harmonic oscillator' )
util.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % Use Fourier grid
space.dof{1}.x_min = -10;                % Lower bound of the grid
space.dof{1}.x_max = 10;                 % Upper bound of the grid
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.mass  = 1;                  % Mass for the kinetic energy

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 020;                   % Index of final time step

time.main.delta  = pi/10;                % Size of time steps 
time.sub.n       =  100;                 % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    =  0.0;           % Lower truncation of energy
hamilt.truncate.max    = 50.0;           % Upper truncation of energy

hamilt.pot.handle = @pot.taylor;         % Taylor series: Harmonic oscillator
hamilt.pot.v{1,1} = [0; 1];              % Force constant

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width = sqrt(1/2)/2;          % Width 
psi.dof{1}.pos_0 = -5.0;                 % Center in position representation
psi.dof{1}.mom_0 =  0.0;                 % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % Contour plot of the Wigner function
plots.density.range.on    = true;        % manually set plotting ranges
plots.density.range.x_min = -8;
plots.density.range.x_max = 8;
plots.density.range.y_min = -8;
plots.density.range.y_max = 8;

plots.expect.energies.max = 20;          % plotting ranges for energy values
