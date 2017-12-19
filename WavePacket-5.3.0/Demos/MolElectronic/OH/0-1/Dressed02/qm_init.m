% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time

util.disp ( '**********************************************' )
util.disp ( 'Vibrationally state selective electronic      ' )
util.disp ( 'excitation of OH radical                      ' )
util.disp ( '                                              ' )
util.disp ( 'Light pulse optimized for a transition        ' )
util.disp ( 'from electronic ground state X2Pi(v=0)        ' )
util.disp ( 'to electronically excited state A2Sigma+ (v=1)' )
util.disp ( '                                              ' )
util.disp ( 'Using 2 dressed states to demonstrate         ' )
util.disp ( 'predominance of 1-photon transition           ' )
util.disp ( '                                              ' )
util.disp ( 'see M.V.Korolkov, G.K.Paramonov               ' )
util.disp ( 'Phys. Rev. A 57(6), 4998  (1998), Fig. 2a+b   ' )
util.disp ( '**********************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation='adi';
hamilt.coupling.labels ={'X ^2\Pi', 'A ^2\Sigma^+'};

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 1728.539;            % Reduced mass (OH radical)
space.dof{1}.n_pts = 064;                % Number of grid points
space.dof{1}.x_min =  0.5;               % Lower bound of grid 
space.dof{1}.x_max =  5.0;               % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 100;                   % Index of final time step

time.main.delta  = 1/atomic.t.fs;        % Size of time steps: 1 fs 
time.sub.n   =   10;                     % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Electric field as sequence of pulses
time.efield.dressed = true;              % Not using dressed state picture
time.efield.photons = {0 -1};            % Photon dressed states
time.efield.shape   = 'sin^2';           % Shape of envelope
time.efield.polar   = 0.0;               % Polarization angle [rad]
time.efield.delay   = 50/atomic.t.fs;    % Time delay of pulse center: 50 fs
time.efield.fwhm    = 50/atomic.t.fs;    % Full width at half maximum: 50 fs
time.efield.ampli   = 142.44/atomic.F.MV_cm; % Field amplitude: 142.44 MV/cm
time.efield.frequ   = 35270.76/atomic.w.cm_1; % Carrier frequency: 35270.76 cm-1
time.efield.phase   = 0.0;               % Phase

% Hamiltonian operator 
hamilt.truncate.min    = -0.2;           % Lower truncation of energy
hamilt.truncate.max    =  0.5;           % Upper truncation of energy

hamilt.pot.handle      = @pot.interp;    % Potentials: Interpolate table

hamilt.dip.handle      = @dip.interp;    % Dipoles: Interpolate table

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % Negative imaginary potential
hamilt.nip.exp  = 4;                     % Exponent
hamilt.nip.min = +1.0;                   % Beginning of inner grid region
hamilt.nip.max = +4.5;                   % End of inner grid region

% Initial wave function
psi.dof{1}.handle = @wav.gauss;          % Gaussian-shaped wavepacket
psi.dof{1}.width  =  0.1331;             % Width 
psi.dof{1}.pos_0  =  1.8601;             % Center in position representation
psi.dof{1}.mom_0  =  0.0;                % Center in momentum representation

psi.init.coeffs         = [1 0];         % Initially: ground state only
psi.init.representation = 'dia';

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % make a contour plot of the Wigner transform

plots.density.range.on = true;           % manually set the ranges for the plot
plots.density.range.x_min = 1;
plots.density.range.x_max = 3.5;
plots.density.range.y_min = -20;
plots.density.range.y_max = 20;

plots.expect.energies.max = 0.05;
plots.expect.energies.min = -0.18;
