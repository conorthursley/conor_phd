% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time

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
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = 0.7;                % Lower bound of grid 
space.dof{1}.x_max = 5.0;                % Upper bound of grid     

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

% Propagator
time.propa.handle = @ket.splitting;      % Operator splitting
time.propa.order = 3;                    % Symmetrized (Strang) scheme

% Electric field loaded from file
time.efield.shape   = 'inter';           % Interpolate from file; We interpolate the
time.efield.ampli   = 1;                 % whole pulse from file, including the
time.efield.delay   = 0;                 % oscillations, so all other variables only get
time.efield.frequ   = 0;                 % dummy values.
time.efield.phase   = 0;
time.efield.polar   = 0;
time.efield.file    = {'tdse_optimal_30.dat'};   % file to load the field from
time.efield.method  = 'spline';          % Spline interpolation

% Initial wave function
psi.dof{1}.handle = @wav.morse;          % Ground state of Morse oscillator
psi.dof{1}.d_e = hamilt.pot.d_e;         % data copied from Morse potential
psi.dof{1}.r_e = hamilt.pot.r_e;         % data copied from Morse potential
psi.dof{1}.alf = hamilt.pot.alf;         % data copied from Morse potential
psi.dof{1}.n_q = 0;                      % ground state

% Modify settings for appearance of plots (if desired)
plots.density.type = 'curve';            % Simple curve plot of 
plots.expect.energies.max = 0.2;         % Range of energy plot
