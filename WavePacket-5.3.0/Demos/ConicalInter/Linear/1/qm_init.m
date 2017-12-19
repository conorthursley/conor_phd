% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '***************************************' )
util.disp ( 'Linear J-T effect: Conical intersection' )
util.disp ( 'diabatic representation                ' )
util.disp ( '***************************************' )

% Number of (coupled) Schrödinger equations
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'dia';

% Spatial discretization
space.dof{1}       = grid.fft;           % using FFT grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = -10.0;              % Lower bound of grid 
space.dof{1}.x_max =  10.0;              % Upper bound of grid

% Spatial discretization
space.dof{2}       = grid.fft;           % using FFT grid
space.dof{2}.mass  = 1;                  % mass
space.dof{2}.n_pts = 64;                 % Number of grid points
space.dof{2}.x_min = -10.0;              % Lower bound of grid 
space.dof{2}.x_max =  10.0;              % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 019;                   % Index of final time step

time.main.delta = 0.10;                  % Size of time steps 
time.sub.n      = 0050;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = -10.0;          % Lower truncation of energy
hamilt.truncate.max    =  50.0;          % Upper truncation of energy

hamilt.pot.handle = @pot.con_int;        % Conical intersection
hamilt.pot.omega=  0.0;                  % Harmonic frequency
hamilt.pot.kappa=  1.0;                  % Linear JT coupling
hamilt.pot.gamma=  0.0;                  % Quadratic JT coupling

hamilt.nip.handle = @nip.power;          % Negative imaginary potential
hamilt.nip.exp  = [4 4];                 % Exponent
hamilt.nip.min  = [-08.0 -08.0];         % Lower bound
hamilt.nip.max  = [+08.0 +08.0];         % Upper bound

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.7;                 % Width 
psi.dof{1}.pos_0 = -4.0;                 % Center in position representation
psi.dof{1}.mom_0 =  4.0;                 % Center in momentum representation

psi.dof{2}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{2}.width =  1.0;                 % Width 
psi.dof{2}.pos_0 =  0.0;                 % Center in position representation
psi.dof{2}.mom_0 =  0.0;                 % Center in momentum representation

psi.init.representation = 'dia';         % Initial WF in diabatic representation
psi.init.coeffs       = [0 1];           % Initially only repulsive diabatic state populated


% Modify settings for appearance of plots (if desired)
plots.density.type = 'surface';

plots.density.surface.view  = [25 25];   % View point for surface plot: [az el]

% The plot ranges look ugly for surface plots
if strcmp(plots.density.type, 'contour')
    plots.density.range.on = true;       % manual setting of plotting range
    plots.density.range.x_min = -7;
    plots.density.range.x_max = 8.5;
    plots.density.range.y_min = -4;
    plots.density.range.y_max = 4;
end

plots.expect.energies.max = 20;          % manual setting of upper energy bond
