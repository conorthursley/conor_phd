% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '******************************************' )
util.disp ( 'Generic E x e conical intersection with   ' )
util.disp ( 'linear and quadratic Jahn-Teller coupling:' )
util.disp ( 'Dynamics on warped, lower state only' )
util.disp ( 'Localized wavepacket with k=5' )
util.disp ( '******************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 064;                % Number of grid points
space.dof{1}.x_min = -7.5;               % Lower bound of grid 
space.dof{1}.x_max =  7.5;               % Upper bound of grid

% Spatial discretization
space.dof{2} = grid.fft;                 % using FFT grid
space.dof{2}.mass = 1;                   % mass
space.dof{2}.n_pts = 064;                % Number of grid points
space.dof{2}.x_min = -7.5;               % Lower bound of grid 
space.dof{2}.x_max =  7.5;               % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 050;                   % Index of final time step

time.main.delta = 0.05;                  % Size of time steps 
time.sub.n      = 0050;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min = -15.0;             % Lower truncation of dia. energy
hamilt.truncate.max = +25.0;             % Upper truncation of dia. energy

hamilt.pot.handle = @pot.con_int;        % Conical intersection
hamilt.pot.omega=  5.0;                  % Harmonic frequency
hamilt.pot.kappa= 10.0;                  % Linear JT coupling
hamilt.pot.gamma=  1.0;                  % Quadratic JT coupling

% Initial wave function
a = 0.60*pi;

psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.5;                 % Width 
psi.dof{1}.pos_0 =  2.5;                 % Center in position representation
psi.dof{1}.mom_0 =  5 * cos(a);          % Center in momentum representation

psi.dof{2}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{2}.width =  0.5;                 % Width 
psi.dof{2}.pos_0 =  0.0;                 % Center in position representation
psi.dof{2}.mom_0 =  5 * sin(a);          % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.density.surface.view  = [45 75];   % Perspective view 

if strcmp(plots.density.type, 'contour')
    plots.density.range.on    = true;    % Manual setting of plotting range
    plots.density.range.x_min = -5;
    plots.density.range.x_max =  5;
    plots.density.range.y_min = -5;
    plots.density.range.y_max =  5;
end

plots.expect.energies.max = 15;
