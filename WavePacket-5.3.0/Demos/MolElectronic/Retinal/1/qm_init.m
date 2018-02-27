% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time

util.disp ( '****************************************************' )
util.disp ( 'Retinal isomerization: model in 1 dimension' )
util.disp ( 'S.Hahn, G.Stock, J. Phys. Chem. B 104(6),1149 (2000)' )
util.disp ( '****************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'adi';

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 56198.3470;          % mass
space.dof{1}.n_pts = 384;                % Number of grid points
space.dof{1}.x_min = -pi;                % Lower bound of grid 
space.dof{1}.x_max =  pi;                % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 050;                   % Index of final time step

time.main.delta  = 100;                  % Size of time steps 
time.sub.n       = 100;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Operator splitting
time.propa.order = 3;                    % Strang method

% Hamiltonian operator 
hamilt.truncate.min    =  -0.1;          % Lower truncation of energy
hamilt.truncate.max    =  +0.5;          % Upper truncation of energy

hamilt.pot.handle = @pot.retinal_1D;     % Hahn&Stock: Retinal potential
hamilt.pot.shift   = [0.0  2.48] / atomic.E.eV; % Vertical offset
hamilt.pot.barrier = [3.6 -1.09] / atomic.E.eV; % Rotational barrier
hamilt.pot.omega   = [0.19 0.19] / atomic.E.eV; % Vibrational frequency
hamilt.pot.kappa   = [0.0, 0.10] / atomic.E.eV; % Linear parameter
hamilt.pot.xc      = 0.5;                       % Slice: constant coupling coordinate
hamilt.pot.lambda  = 0.19 / atomic.E.eV;        % Linear coupling parameter

% Initial wave function
psi.dof{1}.handle = @wav.gauss;          % Gaussian-shaped wavepacket
psi.dof{1}.width  =  0.25;               % Width 
psi.dof{1}.pos_0  =  0;                  % Center in position representation
psi.dof{1}.mom_0  =  0;                  % Center in momentum representation

psi.init.representation = 'adi';
psi.init.coeffs       = [0 1];           % Initially only upper adiabatic state populated

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'contour';   % Contour plot of Wigner transform

plots.expect.energies.max = 0.1;         % Range for energy plot
