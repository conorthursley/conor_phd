% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%                    Boris Schaefer-Bung
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time control reduce

util.disp ( '**********************************' )
util.disp ( 'Asymmetric double well potential  ' )
util.disp ( ' see example 1 in:                ' )
util.disp ( ' B. Schaefer-Bung, C. Hartmann,   ' )
util.disp ( ' B. Schmidt, and Ch. Schuette,    ' )
util.disp ( ' J. Chem. Phys. 135, 014112 (2011)' )
util.disp ( '**********************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 162;                 % I=162 [hbar^2/D], see section IV.B
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -2.0;               % Lower bound of grid 
space.dof{1}.x_max =  2.0;               % Upper bound of grid

% Temporal discretization
time.main.start = 0;                     % Index of initial time step
time.main.stop  = 100;                   % Index of final time step
time.main.delta = 3*pi/(100*0.196);      % Size of time steps 

% Propagator
time.propa.handle = @ket.splitting;      % Operator splitting
time.propa.order = 3;                    % Symmetrized (Strang) scheme

% Electric field as half-cycle pulse
time.efield.dressed = false;             % Not using dressed state picture
time.efield.shape   = 'recta';           % Rectangular shape of envelope
time.efield.polar   = 0;                 % Polarization angle [rad]
time.efield.delay   = pi/(2*0.196);      % Time delay of pulse center
time.efield.fwhm    = pi/0.196;          % Full width at half maximum
time.efield.ampli   = 7.5;               % u_0=7.5 [D/(e*a_0)] see section IV.B
time.efield.frequ   = 0.196;             % Omega = 0.196 [D/hbar], see section IV.B
time.efield.phase   = 0;                 % Phase

% Hamiltonian operator 
hamilt.truncate.min    =  -1.0;          % Lower truncation of energy
hamilt.truncate.max    =   10.0;         % Upper truncation of energy

hamilt.pot.handle = @pot.taylor;         % Taylor series: Double well potential
hamilt.pot.v{1,1} = [0.055;-4;0;+24];    % Linear, quadratic, quartic constant
hamilt.pot.c{1,1} = 1;                   % Constant energy shift

hamilt.dip.handle = @dip.taylor;         % Dipole moment: Taylor series
hamilt.dip.x.v{1,1} = 1;                 % Linear dipole moment, slope 1

hamilt.sbc.handle = @sbc.taylor;         % System-bath coupling: Taylor series
hamilt.sbc.v{1,1} = 1;                   % Linear coupling, slope 1

% Calculate and save (bound) eigen states (==> qm_bound)
psi.eigen.start        = 0;              % Lower index
psi.eigen.stop         = 20;             % Upper index
psi.save.export        = true;           % Toggle saving
psi.save.file = 'bound';                 % Filename for saving
psi.save.dir = pwd;                      % Directory for saving

% Select ODE solver, parameters (==> qm_control)
control.solvers.handle2 = @ode45;        % Runge-Kutta
control.solvers.reltol = 1e-6;           % Relative error tolerance for ode45
          
% Initial density
control.initial.choice = 'thermal';      % thermal = Boltzmann density
control.initial.temperature = 1.0;       % choice of temperature

% Temperature and system-bath coupling
control.lvne.order = 'df';               % ordering vectorized density matrices: diagonals first
control.lvne.temperature = 1.0;          % temperature in atomic units: 315,923.5 K
control.relax.rate =  0.264386355;       % Gamma_{2->0} should be = 1 [D/hbar]
                                         % set to 0.2644 due to error in ancestor code 
control.relax.lower = 0;                 % Lower state for reference transition
control.relax.upper = 2;                 % Upper state for reference transition
control.relax.model = 'fermi';           % Omega dependent (Andrianov&Saalfrank)

% Define output observables and choose control targets
control.observe.types = 'prj';           % choose populations of observables
control.observe.choices = {0:2:10 1:2:9 11:20}; 
control.observe.labels = {'left well' 'right well' 'delocalized'};
control.observe.targets = 1:3;           % choose control targets 

% Parameters for balancing transformation
reduce.balance.A_stable = 'ssu';         % SSU or EVS method for stabilizing A  
reduce.balance.A_shift = 1e-6;           % magnitude of eigenvalue shift
reduce.balance.A_split = 1;              % dimension of instable part of A
reduce.balance.BN_scale = 3;             % Scaling factor for control fields
reduce.balance.acf_couple = true;        % additional control field coupling A and rho/x
reduce.balance.method = 'iter';          % ITER or BICG solver for GLYAP
reduce.balance.transform = 'srbt';       % SRBT or MRMR balancing transform;

% Parameters for H2 optimal model reduction (Benner/Breiten @ MPI Magdeburg)
reduce.H2model.A_stable = 'ssu';         % SSU or EVS method for stabilizing A  
reduce.H2model.A_shift = 1e-6;           % magnitude of eigenvalue shift
reduce.H2model.A_split = 1;              % dimension of instable part of A
reduce.H2model.BN_scale = 3;             % Scaling factor for control fields
reduce.H2model.acf_couple = true;        % additional control field coupling A and rho/x
reduce.H2model.method = 'iter';          % ITER or BICG solver for BIRKA

% Parameters for calculation of H2 error
reduce.H2error.A_shift = 1e-2;           % magnitude of eigenvalue shift
reduce.H2error.BN_scale = 10;            % Scaling factor for control fields
reduce.H2error.method = 'iter';          % ITER or BICG solver for GSYLV

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % Contour plot of Wigner transform
plots.density.export.on = true;
plots.expect.export.on = true;
plots.expect.population.min = 0.26;
plots.expect.population.max = 0.42;

plots.expect.energies.max = 2;           % Range of energy plot


