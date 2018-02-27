% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '********************************************' )
util.disp ( 'Spin-boson system in 2 dimensions: adiabatic' )
util.disp ( '********************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'adi';

% Spatial discretization
space.dof{1}       = grid.fft;           % using FFT grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.n_pts = 192;                % Number of grid points
space.dof{1}.x_min = -12.0;              % Lower bound of grid 
space.dof{1}.x_max =  12.0;              % Upper bound of grid

% Spatial discretization
space.dof{2}       = grid.fft;           % using FFT grid
space.dof{2}.mass  = 1;                  % mass
space.dof{2}.n_pts = 192;                % Number of grid points
space.dof{2}.x_min = -12.0;              % Lower bound of grid 
space.dof{2}.x_max =  12.0;              % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 042;                   % Index of final time step

time.main.delta = 0.1;                   % Size of time steps 
time.sub.n      = 0150;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.cheby_real;     % Expansion in Chebychev polynomials
time.propa.order = 0;                    % auto-determine order
time.propa.precision = 1e-10;            % required precision

% Hamiltonian operator 
hamilt.truncate.min    = -500.0;         % Lower truncation of dia. energy
hamilt.truncate.max    = +500.0;         % Upper truncation of dia. energy

hamilt.pot.handle = @pot.con_int;        % Conical intersection
hamilt.pot.omega= 03.0;                  % Harmonic frequency
hamilt.pot.kappa= 06.0;                  % Linear JT coupling
hamilt.pot.gamma= 00.0;                  % Quadratic JT coupling

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.5;                 % Width 
psi.dof{1}.pos_0 = -5.0;                 % Center in position representation
psi.dof{1}.mom_0 =  0.0;                 % Center in momentum representation

psi.dof{2}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{2}.width =  2.5;                 % Width 
psi.dof{2}.pos_0 =  0.0;                 % Center in position representation
psi.dof{2}.mom_0 =  0.0;                 % Center in momentum representation

psi.init.representation = 'dia';
psi.init.coeffs       = [0 1];           % Initially left diabatic state populated

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.density.surface.view = [15 45];   % Perspective view (with potentials)

plots.expect.energies.max = 150;
