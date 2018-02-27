% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '************************************************' )
util.disp ( 'Squeezed/coherent state of a harmonic oscillator' )
util.disp ( '************************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % Using FFT grid 
space.dof{1}.n_pts = 032;                % Number of grid points
space.dof{1}.x_min = 05.0;               % Lower bound of grid 
space.dof{1}.x_max = 15.0;               % Upper bound of grid
space.dof{1}.mass = 1.0;                 % Mass

space.dof{2} = grid.fft;                 % Using FFT grid 
space.dof{2}.n_pts = 064;                % Number of grid points
space.dof{2}.x_min = 15.0;               % Lower bound of grid 
space.dof{2}.x_max = 25.0;               % Upper bound of grid
space.dof{2}.mass = 1.0;                 % Mass

space.dof{3} = grid.fft;                 % Using FFT grid 
space.dof{3}.n_pts = 064;                % Number of grid points
space.dof{3}.x_min = 25.0;               % Lower bound of grid 
space.dof{3}.x_max = 35.0;               % Upper bound of grid
space.dof{3}.mass = 1.0;                 % Mass

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 050;                   % Index of final time step

time.main.delta = pi/25;                 % Size of time steps 
time.sub.n      = 0100;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    =   0.0;          % Lower truncation of energy
hamilt.truncate.max    = 100.0;          % Upper truncation of energy

hamilt.pot.handle      = @pot.taylor;    % Taylor series: Harmonic oscillator
hamilt.pot.v{1,1} = [0 0 0;1 1 1];       % Force constants
hamilt.pot.s{1,1} = [10  20  30 ];       % Minimum positions


% Periodic boundary conditions
hamilt.nip.handle = @nip.power;          % Negativ imaginary potential
hamilt.nip.exp  = [2  2  2 ];            % Exponent
hamilt.nip.min  = [02 12 22];            % Lower bound
hamilt.nip.max  = [18 28 38];            % Upper bound

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width = 1.5*sqrt(1/2);        % Width: squeezed state
psi.dof{1}.pos_0 = 10.0;                 % Center in position representation
psi.dof{1}.mom_0 =  0.0;                 % Center in momentum representation

psi.dof{2}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{2}.width = sqrt(1/2);            % Width: coherent state
psi.dof{2}.pos_0 =  22.5;                % Center in position representation
psi.dof{2}.mom_0 =  +0.0;                % Center in momentum representation

psi.dof{3}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{3}.width = sqrt(1/2);            % Width: coherent state
psi.dof{3}.pos_0 =  30.0;                % Center in position representation
psi.dof{3}.mom_0 =  -2.5;                % Center in momentum representation


% Modify settings for appearance of plots (if desired)
plots.density.type          = 'surface'; % 3D plot
plots.density.surface.view  = [60 45];   % View point for surface plot: [az el]
plots.density.representation = 'dvr';    % DVR or FBR representation
plots.expect.energies.max = 10;          % Manually set range for energies
