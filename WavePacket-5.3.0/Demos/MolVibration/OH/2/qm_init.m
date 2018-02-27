% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic control hamilt psi space time

util.disp ( '****************************************************' )
util.disp ( 'Vibrationally state selective excitation of a Morse ' )
util.disp ( 'oscillator resembling one of the OH bonds in water  ' )
util.disp ( 'Laser pulse optimized for 5-photon transition: 0->5 ' )
util.disp ( '                                                    ' )
util.disp ( 'see M.V.Korolkov, G.K.Paramonov, and B. Schmidt     ' )
util.disp ( 'J. Chem. Phys. 105(5), 1862-1879  (1996), Fig. 2a   ' )
util.disp ( '****************************************************' )

% In the paper we used a reduced mass of 1728.539, but here 1728.257
m_O = atomic.mass.O16;
m_H = atomic.mass.H1;

% Spatial discretization
space.dof{1} = grid.fft;                 % Using fft grid
space.dof{1}.mass  = m_O*m_H/(m_O+m_H);  % Reduced mass 
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min = 0.7;                % Lower bound of grid 
space.dof{1}.x_max = 12.0;               % Upper bound of grid

% Hamiltonian operator 
hamilt.truncate.min = -0.1;              % Lower truncation of energy
hamilt.truncate.max =  0.5;              % Upper truncation of energy

hamilt.pot.handle = @pot.morse;          % Morse oscillator
hamilt.pot.d_e = 0.1994;                 % Dissociation energy
hamilt.pot.r_e = 1.821;                  % Equilibrium length
hamilt.pot.alf = 1.189;                  % Range parameter

hamilt.dip.handle = @dip.mecke;          % Mecke dipole function
hamilt.dip.r_0 = 0.6/atomic.l.A;         % Length parameter: 0.6 A
hamilt.dip.q_0 = 7.85/atomic.d.D*atomic.l.A; % Charge parameter: 7.85 D/A

% Temporal discretization
time.main.delta = 1/atomic.t.fs;         % Size of time steps: 10 fs 
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 1000;                  % Index of final time step
time.sub.n      = 10;                    % Number of sub steps per time step

% Calculate and save (bound) eigen states (==> qm_bound)
psi.eigen.stop = 00;                     % Lower index
psi.eigen.stop = 21;                     % Upper index
psi.save.export = true;                  % Toggle saving
psi.save.file = 'bound';                 % Filename for saving
psi.save.dir = pwd;                      % Directory for saving

% Electric field pulses
time.efield.shape = 'sin^2';             % Shape of envelope
time.efield.polar = 0.0;                 % Polarization angle [rad]
time.efield.delay = 500/atomic.t.fs;     % Time delay of pulse center
time.efield.fwhm  = 500/atomic.t.fs;     % Full width at half maximum
time.efield.ampli = 328.5/atomic.F.MV_cm;% From MV/cm to E_H/(e*a_0)
time.efield.frequ = 3424.19/atomic.w.cm_1;% From cm-1 to E_H
time.efield.phase = 0;                   % Phase

% Initial wave function
% Define populations as observables
control.observe.types = 'prj';           % types of observables
control.observe.choices = {[0] [1] [2] [3] [4] [5] [6] [7:21]};
control.observe.labels  = {'|0>' '|1>' '|2>' '|3>' '|4>' '|5>' '|6>' '|7:21>'};
control.observe.targets = 1:8;           % evaluate targets at terminal time

% Initial state (==> qm_abncd)
control.initial.choice = 'pure';         % starting from a pure state
control.initial.pure = 0;                % vibrational ground state

% Propagator (==> qm_control)
control.solvers.handle2 = @ode45;        % Runge-Kutta: Dormand-Prince (built-in)
control.solvers.reltol = 1e-6;           % Relative error tolerance of solution
