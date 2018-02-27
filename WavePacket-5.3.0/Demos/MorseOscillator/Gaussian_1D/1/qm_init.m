% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '*************************************' )
util.disp ( 'Gaussian packet in a Morse oscillator' )
util.disp ( '*************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1728.539;            % Reduced mass (OH molecule)
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = 0.7;                % Lower bound of grid 
space.dof{1}.x_max = 7.0;                % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 300;                   % Index of final time step

time.main.delta  = 05;                   % Size of time steps 
time.sub.n   =    020;                   % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    =  0.0;           % Lower truncation of energy
hamilt.truncate.max    =  1.0;           % Upper truncation of energy

hamilt.pot.handle      = @pot.morse;     % Harmonic oscillator
hamilt.pot.d_e  = 0.1994;                % Dissociation energy
hamilt.pot.r_e  = 1.821;                 % Equilibrium length
hamilt.pot.alf  = 1.189;                 % Range parameter
hamilt.pot.t_e  = 0.0;                   % Energetic shift

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % Negativ imaginary potential
hamilt.nip.exp  = 4;                     % Exponent
hamilt.nip.min  = 1.0;                   % Beginning of non-gobbler region
hamilt.nip.max  = 6.0;                   % End of non-gobbler region

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.1;                 % Width 
psi.dof{1}.pos_0 =  1.44;                % Center in position representation
psi.dof{1}.mom_0 =  0.0;                 % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'surface';   % Contour plot of Wigner transform

plots.expect.energies.max = 0.1;         % Manually set range for energy plot
