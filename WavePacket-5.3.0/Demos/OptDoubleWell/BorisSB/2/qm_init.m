% Copyright (C) 2004-2017 Burkhard Schmidt's group
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
util.disp ( ' see:                ' )
util.disp ( ' B. Schaefer-Bung, C. Hartmann,   ' )
util.disp ( ' B. Schmidt, and Ch. Schuette,    ' )
util.disp ( ' J. Chem. Phys. 135, 014112 (2011)' )
util.disp ( '**********************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 162;                 % I=162 [hbar^2/D]
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -2.0;               % Lower bound of grid 
space.dof{1}.x_max =  2.0;               % Upper bound of grid

% Temporal discretization
time.main.start = 0;                     % Index of initial time step
time.main.stop  = 050;                   % Index of final time step
time.main.delta = 002;                   % Size of time steps 
time.sub.n      = 020;                   % Number of substeps per main step

% Electric field as half-cycle pulse
time.efield.dressed = false;             % Using "bare state" picture
time.efield.shape   = 'recta';           % Shape of envelope
time.efield.polar   = 0;                 % Polarization angle [rad]
time.efield.delay   = 50;                % Time delay of pulse center
time.efield.fwhm    = 100;               % Full width at half maximum
time.efield.ampli   = 1.00;              % Amplitude
time.efield.frequ   = 0;                 % No carrier frequency
time.efield.phase   = 0;                 % Phase

% frequency resolved optical gating
time.frog.choice    = 'Gauss';           % choice of methods
time.frog.width     = 30;                % width of Gaussian
time.frog.zoom      = 30;                % zooming frequency axis

% Hamiltonian operator 
hamilt.eigen.cutoff    =  0.0;           % Cut-off of Hamiltonian matrix
hamilt.truncate.min    =  -1.0;          % Lower truncation of energy
hamilt.truncate.max    =   10.0;         % Upper truncation of energy

hamilt.pot.handle = @pot.taylor;         % Taylor series: Double well potential
hamilt.pot.v{1,1} = [0.055;-4;0;+24];    % Linear, quadratic, quartic constant
hamilt.pot.c{1,1} = 1;                   % Constant energy shift

hamilt.dip.handle = @dip.taylor;         % Dipole moment: Taylor series
hamilt.dip.x.v{1,1} = 1;                 % Linear dipole moment, slope 1

hamilt.sbc.handle = @sbc.taylor;         % System-bath coupling: Taylor series
hamilt.sbc.v{1,1} = 1;                   % Linear coupling, slope 1

% Calculate and save (bound) eigen states 
psi.eigen.start = 0;                     % Lower index
psi.eigen.stop = 20;                     % Upper index
psi.save.export = true;                  % Export wavefunctions to data file
psi.save.file = 'bound';                 % Filename for saving
psi.save.dir = pwd;                      % Directory for saving

% Select ODE solvers, parameters
control.solvers.handle1 = @ode.RuKu6;    % Runge-Kutta or other integrators
control.solvers.handle2 = @ode45;        % Runge-Kutta: Dormand-Prince (built-in)
control.solvers.reltol = 1e-6;           % Relative error tolerance of solution

% Select initial state (==> qm_control)
control.initial.choice = 'thermal';      % thermal = Boltzmann density
control.initial.temperature = 0.1;       % temperature in atomic units: 315,923.5 K

% Temperature and system-bath coupling (some parameters specified outside)
control.lvne.order = 'df';               % vectorize matrices: diagonals first
control.lvne.temperature = control.initial.temperature;          % ratio of upwards/downwards rates
control.relax.rate =  1/100;             % Gamma_{2->0} 
control.relax.lower = 0;                 % Lower state for reference transition
control.relax.upper = 2;                 % Upper state for reference transition
control.relax.model = 'fermi';           % Omega dependent (Andrianov&Saalfrank)

% Define output observables and choose control targets
control.observe.types = 'prj';           % choose populations of observables
control.observe.choices ={0:2:10 1:2:9 11:20}; 
control.observe.labels ={'left well' 'right well' 'delocalized'};
control.observe.targets = 1:3;           % choose control targets 

% Optimal control theory
control.optimal.terminal = 2;            % Which observable to be optimized
control.optimal.tolerance = 1e-14;       % Threshold terminating iteration
control.optimal.max_iter = 025;          % Max. number of iteration steps
control.optimal.alpha = 0.010;           % Penalty factor for laser fluence
control.optimal.eta  = 1.0;              % Calculate optimal backward field
control.optimal.zeta = 1.0;              % Calculate optimal forward field
control.optimal.order = 2;               % Error order for optimal fields
control.optimal.fb_test = false;         % Only testing forward/backward
control.plot.mov = true;                 % Animation of u(t), x(t), y(t)
control.plot.uxy = true;                 % Plot u(t), x(t), y(t)
control.plot.j12 = true;                 % Plot j_1(t), j_2(t), and total
control.plot.psd = true;                 % Plot power spectral density

% Parameters for balancing transformation
reduce.balance.A_stable = 'ssu';         % SSU or EVS method for stabilizing A 
reduce.balance.A_split = 1;              % dimension of instable part of A
reduce.balance.A_shift = 1e-4;           % value for eigenvalue shift
reduce.balance.BN_scale = 6;             % Scaling factor: divide B,N matrices by this
reduce.balance.acf_couple = true;        % coupling terms between A and rho
reduce.balance.method = 'iter';          % ITER or BICG solver for GLYAP
reduce.balance.transform = 'srbt';       % SRBT or MRMR balancing transform

% Parameters for H2 model reduction (Benner/Breiten @ MPI Magdeburg)
reduce.H2model.A_stable = 'ssu';         % SSU or EVS method for stabilizing A  
reduce.H2model.A_shift = 1e-4;           % magnitude of eigenvalue shift
reduce.H2model.A_split = 1;              % dimension of instable part of A
reduce.H2model.BN_scale = 3;             % Scaling factor for control fields
reduce.H2model.acf_couple = true;        % additional control field coupling A and rho/x
reduce.H2model.method = 'iter';          % ITER or BICG solver for GSYLV
reduce.H2model.max_iter = 200;           % Max. number of iterations
reduce.H2model.conv_tol = 1e-6;          % Convergence tolerance

% Parameters for calculation of H2 error
reduce.H2error.A_shift = 1e-2;           % magnitude of eigenvalue shift
reduce.H2error.BN_scale = 3;             % Scaling factor for control fields
reduce.H2error.method = 'iter';          % ITER or BICG solver for GLYAP

% Appeaarance of various plots
plots.expect.on = false;
plots.expect.energies.max = 4;           % Range of energy plot


