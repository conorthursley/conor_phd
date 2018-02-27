% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '***************************************************' )
util.disp ( 'Dual crossing example, k=30 (Tully 1990)' )
util.disp ( '***************************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'adi';

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 2000;                % mass
space.dof{1}.n_pts = 768;                % Number of grid points
space.dof{1}.x_min = -12.0;              % Lower bound of grid 
space.dof{1}.x_max =  12.0;              % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 040;                   % Index of final time step

time.main.delta  = 025.0;                % Size of time steps 
time.sub.n       = 050;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min = -0.05;             % Lower truncation of energy
hamilt.truncate.max =  1.15;             % Upper truncation of energy

hamilt.pot.handle = @pot.tully2;         % Double Crossing
hamilt.pot.A = 0.10;                     % Tully's parameter A
hamilt.pot.B = 0.28;                     % Tully's parameter B
hamilt.pot.C = 0.015;                    % Tully's parameter C
hamilt.pot.D = 0.06;                     % Tully's parameter D
hamilt.pot.E = 0.05;                     % Tully's parameter E

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % Negative imaginary potential
hamilt.nip.exp  = 4;                     % Exponent
hamilt.nip.min = -10;                    % Beginning of inner grid region
hamilt.nip.max = +10;                    % End of inner grid region

% Initial wave function
psi.dof{1}.handle = @wav.gauss;          % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.75;                % Width 
psi.dof{1}.pos_0   = -7.0;               % Center in position representation
psi.dof{1}.mom_0   = 30.0;               % Center in momentum representation

psi.init.representation = 'adi';
psi.init.coeffs       = [1 0];           % Initially only lower adiabatic state 

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.density.surface.view  = [80 06];   % View point for surface plot: [az el]

if strcmp(plots.density.type, 'contour')
    plots.density.range.on = true;       % manual setting of plotting range
    plots.density.range.x_min = -9;
    plots.density.range.x_max = 9;
    plots.density.range.y_min = 10;
    plots.density.range.y_max = 40;
end

plots.expect.energies.max = 0.3;

