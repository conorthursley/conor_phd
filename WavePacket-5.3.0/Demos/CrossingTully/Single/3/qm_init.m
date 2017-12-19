% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '***************************************************' )
util.disp ( 'Single crossing example, k=10 (Tully 1990)' )
util.disp ( '***************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'adi';

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 2000;                % mass
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  10.0;              % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 020;                   % Index of final time step

time.main.delta  = 100.0;                % Size of time steps 
time.sub.n       = 050;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = -0.01;          % Lower truncation of energy
hamilt.truncate.max    =  0.10;          % Upper truncation of energy

hamilt.pot.handle = @pot.tully1;         % Single Crossing
hamilt.pot.A = 0.01;                     % Tully's parameter A
hamilt.pot.B = 1.6;                      % Tully's parameter B
hamilt.pot.C = 0.005;                    % Tully's parameter C
hamilt.pot.D = 1.0;                      % Tully's parameter D

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % Negative imaginary potential
hamilt.nip.exp  = 4;                     % Exponent
hamilt.nip.min = -9;                     % Beginning of inner grid region
hamilt.nip.max =  9;                     % End of inner grid region

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.75;                % Width 
psi.dof{1}.pos_0 = -6.0;                 % Center in position representation
psi.dof{1}.mom_0 = 10.0;                 % Center in momentum representation

psi.init.representation = 'adi';
psi.init.coeffs       = [1 0];           % Initially only lower adiabatic state 

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.density.surface.view  = [84 06];   % View point for surface plot: [az el]

if strcmp(plots.density.type, 'contour')
    plots.density.range.on = true;       % manual setting of plotting range
    plots.density.range.x_min = -10;
    plots.density.range.x_max = 10;
    plots.density.range.y_min = 0;
    plots.density.range.y_max = 15;
end

plots.expect.energies.max = 0.04;

% The upper wave function is cut off by default, so tune the displayed energy range
plots.density.pot.max = 0.015;
