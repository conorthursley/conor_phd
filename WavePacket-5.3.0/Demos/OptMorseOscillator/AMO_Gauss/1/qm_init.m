% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time control 

util.disp ( '****************************************************' )
util.disp ( 'Optimal control of a bond length for a model Morse  ' )
util.disp ( 'oscillator resembling one of the OH bonds in water  ' )
util.disp ( '                                                    ' )
util.disp ( 'see Figs. 1-4 of: W. Zhu and H. Rabitz              ' )
util.disp ( 'J. Chem. Phys. 109(2), 385-391  (1998)              ' )
util.disp ( 'http://dx.doi.org/10.1063/1.476575                  ' )
util.disp ( '****************************************************' )

% Isotopic masses
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization ==> qm_bound
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass (OH radical)
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min =  0.7;               % Lower bound of grid 
space.dof{1}.x_max = 10.0;               % Upper bound of grid     

space.amo{1}.handle = @amo.gauss;        % Choose Gaussian-shaped AMO
space.amo{1}.pos_0 = 2.5;                % Center of Gaussian
space.amo{1}.width = 1/(2*25);           % Width parameter of Gaussian
space.amo{1}.label = 'Gaussian|2.5|25';  % Labeling this AMO

% Hamiltonian operator 
hamilt.truncate.min = -0.1;              % Lower truncation of energy
hamilt.truncate.max = +1.0;              % Upper truncation of energy

hamilt.pot.handle = @pot.morse;          % Morse oscillator
hamilt.pot.d_e = 0.1994;                 % Dissociation energy
hamilt.pot.r_e = 1.821;                  % Equilibrium length
hamilt.pot.alf = 1.189;                  % Range parameter

hamilt.dip.handle = @dip.mecke;          % Mecke dipole function
hamilt.dip.r_0 = 0.6/atomic.l.A;         % Length parameter: 0.6 A
hamilt.dip.q_0 = 7.85/atomic.d.D*atomic.l.A; % Charge parameter: 7.85 D/A

% Temporal discretization
time.main.delta = 10/atomic.t.fs;        % Size of time steps: 10 fs 
time.main.start = 0000;                  % Index of initial time step
time.main.stop  = 0317;                  % Index of final time step
time.sub.n      = 200;                   % Number of sub steps per time step

% Calculate and save (bound) eigen states (==> qm_bound)
psi.eigen.stop = 00;                     % Lower index
psi.eigen.stop = 21;                     % Upper index
psi.save.export = true;                  % Toggle saving
psi.save.file = 'bound';                 % Filename for saving
psi.save.dir = pwd;                      % Directory for saving

% Electric field pulses: Initial guess
time.efield.shape   = 'recta';           % Shape of envelope
time.efield.polar   = 0;                 % Polarization angle [rad]
time.efield.delay   = 1585/atomic.t.fs;  % Time delay of pulse center
time.efield.fwhm    = 3170/atomic.t.fs;  % Pulse length
time.efield.ampli   = 0.02;              % Amplitude of electric field         
time.efield.frequ   = 0;                 % Carrier frequency
time.efield.phase   = 0;                 % Phase

time.frog.choice = 'none';               % Frequency resolved optical gating
time.frog.zoom = 15;                     % Zoom factor for frequency axes 

% Select ODE solver, parameters
control.solvers.handle1 = @ode.RuKu4;    % Runge-Kutta or other integrators
control.solvers.handle2 = @ode45;        % Runge-Kutta: Dormand-Prince (built-in)
control.solvers.reltol = 1e-6;           % Relative error tolerance of solution

% Initial state (==> qm_abncd)
control.initial.choice = 'pure';         % starting from a pure state
control.initial.pure = 0;                % vibrational ground state

% Define outputs: projections as observables
control.observe.types = 'amo';           % types of observables
control.observe.choices = {1};
control.observe.targets = 1;             % choose control targets 

% Optimal control 
control.optimal.terminal = 1;            % Which observable to be optimized
control.optimal.tolerance = 1e-20;       % Threshold terminating iteration
control.optimal.max_iter = 030;          % Max. number of iteration steps
control.optimal.alpha = 1.00;            % Penalty factor for laser fluence
control.optimal.eta  = 1.00;             % Calculate optimal backward field
control.optimal.zeta = 1.00;             % Calculate optimal forward field
control.optimal.order = 2;               % Error order for optimal fields
control.optimal.prefactor = 'current';   % How to calculate prefactor (C2)
control.optimal.fb_test = false;         % Only testing forward/backward

% Plotting options
control.plot.uxy = true;                 % Plot u(t), x(t), y(t)
control.plot.j12 = true;                 % Plot j_1(t), j_2(t), and total
control.plot.psd = true;                 % Plot power spectral density
control.plot.mov = true;                 % Animation of u(t), x(t), y(t)

% Modify settings for appearance of plots (if desired)
plots.density.type = 'curve';            % Simple curve plot of 
plots.expect.energies.max = 0.2;         % Range of energy plot
