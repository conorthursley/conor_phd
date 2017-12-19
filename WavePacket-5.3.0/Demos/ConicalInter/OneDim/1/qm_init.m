% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '********************************************' )
util.disp ( 'Spin-boson system in 1 dimensions' )
util.disp ( 'Parameters from Monika Hejjas thesis' )
util.disp ( 'HU Berlin, Dept. of Physics, 2004' )
util.disp ( '********************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'adi';

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 1024;               % Number of grid points
space.dof{1}.x_min = -90.0;              % Lower bound of grid 
space.dof{1}.x_max =  120.0;              % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 150;                   % Index of final time step

time.main.delta = 0.5;                   % Size of time steps 
time.sub.n      = 0200;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Splitting method
time.propa.order = 3;                    % Strang

% Hamiltonian operator 
hamilt.truncate.min    = -10.0;          % Lower truncation of dia. energy
hamilt.truncate.max    = +50.0;          % Upper truncation of dia. energy

hamilt.pot.handle = @pot.con_int;        % for conical intersections
hamilt.pot.omega = 00.01;                % Harmonic frequency
hamilt.pot.gamma = 00.50;                % spin coupling
hamilt.pot.kappa = sqrt(0.05);           % asymmetry parameter
hamilt.pot.delta = 2*7.44;               % energy gap

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width = 2.236068;             % Width 
psi.dof{1}.pos_0 = -71.5;                % Center in position representation
psi.dof{1}.mom_0 =  00.0;                % Center in momentum representation

psi.init.representation = 'dia';
psi.init.coeffs       = [1 0];           % Initially left diabatic state populated

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.expect.energies.max = 20;
