% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time

util.disp ( '***************************************' )
util.disp ( 'Two state / three mode model of        ' )
util.disp ( 'S2->S1 internal conversion of pyrazine ' )
util.disp ( '***************************************' )

omega = [0.126 0.074 0.118]/atomic.E.eV; % Harmonic frequencies in eV

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'dia';

% Grid definition
space.dof{1}       = grid.fft;           % Using a FFT grid
space.dof{1}.x_min = -10;                % Lower bound of grid
space.dof{1}.x_max = 10;                 % Upper bound of grid
space.dof{1}.n_pts = 96;                 % Number of grid points
space.dof{1}.mass  = 1/omega(1);         % "Mass" for kinetic energy

space.dof{2}       = grid.fft;           % Using a FFT grid
space.dof{2}.x_min = -10;                % Lower bound of grid
space.dof{2}.x_max = 10;                 % Upper bound of grid
space.dof{2}.n_pts = 96;                 % Number of grid points
space.dof{2}.mass  = 1/omega(2);         % "Mass" for kinetic energy

space.dof{3}       = grid.fft;           % Using a FFT grid
space.dof{3}.x_min = -10;                % Lower bound of grid
space.dof{3}.x_max = 10;                 % Upper bound of grid
space.dof{3}.n_pts = 96;                 % Number of grid points
space.dof{3}.mass  = 1/omega(3);         % "Mass" for kinetic energy

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 300;                   % Index of final time step

time.main.delta = 1/atomic.t.fs;         % Size of time steps: 1 femtosecond 
time.sub.n      =  100;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    =    0.0;         % Lower truncation of energy
hamilt.truncate.max    =   +1.0;         % Upper truncation of energy

hamilt.pot.handle = @pot.pyrazine_3d;    % Pyrazine model potential
hamilt.pot.omega  = omega;               % Force constants
hamilt.pot.kappa1 = [+0.037 -0.105]/atomic.E.eV;  % Linear parameter (tuning)
hamilt.pot.kappa2 = [-0.254 +0.149]/atomic.E.eV;  % Linear parameter (tuning)
hamilt.pot.lambda =  +0.262/atomic.E.eV; % Interstate coupling
hamilt.pot.energ1 =   3.940/atomic.E.eV; % Vertical excitation energy
hamilt.pot.energ2 =   4.840/atomic.E.eV; % Vertical excitation energy

% Initial wave function
psi.corr.handle = @wav.gauss;            % Gaussian-shaped wavepacket
psi.corr.width  = [1 1 1]/sqrt(2);       % Width of coherent state
psi.corr.pos_0  = [0.0 0.0 0.0];         % Center in position representation
psi.corr.mom_0  = [0.0 0.0 0.0];         % Center in momentum representation

psi.init.representation = 'dia';       
psi.init.coeffs         = [0 1];         % Initially upper adiabatic state populated

% Modify settings for appearance of plots (if desired)
plots.density.type  = 'reduced';         % View point for surface plot: [az el]
plots.density.contour.nlev  = [60 15];   % Number of contours: density/energy
