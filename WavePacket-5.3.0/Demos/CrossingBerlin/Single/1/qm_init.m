% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '***************************************************' )
util.disp ( 'Single crossing example (Nettsheim, Schuette: 1999)' )
util.disp ( '***************************************************' )

% Number of (coupled) Schrödinger equations and representation
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'adi';

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 10000;               % mass
space.dof{1}.n_pts = 512;                % Number of grid points
space.dof{1}.x_min = 0.05;               % Lower bound of grid 
space.dof{1}.x_max =  2.85;              % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 100;                   % Index of final time step

time.main.delta  = 1.0;                  % Size of time steps 
time.sub.n       = 050;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min =  0.0;              % Lower truncation of energy
hamilt.truncate.max = 15.0;              % Upper truncation of energy

hamilt.pot.handle = @pot.single;         % Single Crossing
hamilt.pot.repuls  = 1.0;                % Repulsion (exponential)
hamilt.pot.attrac  = 1.0;                % Attraction (parabolic)
hamilt.pot.couple  = 0.1;                % Coupling (constant)

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % Negative imaginary potential
hamilt.nip.exp  = 4;                     % Exponent
hamilt.nip.min = +0.15;                  % Beginning of inner grid region
hamilt.nip.max = +2.50;                  % End of inner grid region

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.05;                % Width 
psi.dof{1}.pos_0 =  0.4;                 % Center in position representation
psi.dof{1}.mom_0 =  100;                 % Center in momentum representation

psi.init.representation = 'adi';
psi.init.coeffs       = [0 1];           % Initially only upper adiabatic state populated

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.density.surface.view  = [70 35];   % View point for surface plot: [az el]

if strcmp(plots.density.type, 'contour')
    plots.density.range.on    = true;    % manual setting of plotting range
    plots.density.range.x_min = 0.2;
    plots.density.range.x_max = 2.7;
    plots.density.range.y_min = -50;
    plots.density.range.y_max = 280;
end

plots.expect.energies.max = 5;
