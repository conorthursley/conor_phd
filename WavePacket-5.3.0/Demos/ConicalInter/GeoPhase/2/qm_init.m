% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '******************************************' )
util.disp ( 'Generic E x e conical intersection example' )
util.disp ( 'with linear Jahn-Teller coupling included:' )
util.disp ( 'Dynamics on lower (Mexican hat) state only' )
util.disp ( 'Neglect of geometric (Berry) phase effect!' )
util.disp ( '******************************************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.x_min = -4.5;               % Lower bound of grid 
space.dof{1}.x_max =  4.5;               % Upper bound of grid

% Spatial discretization
space.dof{2} = grid.fft;                 % using fft grid
space.dof{2}.mass = 1;                   % mass
space.dof{2}.n_pts = 64;                 % Number of grid points
space.dof{2}.x_min = -4.5;               % Lower bound of grid 
space.dof{2}.x_max =  4.5;               % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 100;                   % Index of final time step

time.main.delta = 0.05;                  % Size of time steps 
time.sub.n      = 0050;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Splitting method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = -200.0;         % Lower truncation of dia. energy
hamilt.truncate.max    = +200.0;         % Upper truncation of dia. energy

hamilt.pot.handle = @pot.con_int;        % Conical intersection
hamilt.pot.omega=  5.0;                  % Harmonic frequency
hamilt.pot.kappa= 10.0;                  % Linear JT coupling
hamilt.pot.gamma=  0.0;                  % Quadratic JT coupling

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.5;                 % Width 
psi.dof{1}.pos_0 = -2.0;                 % Center in position representation
psi.dof{1}.mom_0 =  0.0;                 % Center in momentum representation

psi.dof{2}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{2}.width =  0.5;                 % Width 
psi.dof{2}.pos_0 =  0.0;                 % Center in position representation
psi.dof{2}.mom_0 =  0.0;                 % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type = 'surface';

plots.density.surface.view  = [30 60];   % Side view 
plots.density.energy.on     = true;      % With energy surfaces, densities color-coded

% plots.density.surface.view  = [90 -90];  % Bottom view
% plots.density.energy.on     = false;     % No energy surfaces, only densities
% plots.density.surface.look  = [ 1  1];   % Neither shading nor lighting
% plots.density.surface.lite  = [90 -90];  % Angle of light for surface plot: [az el]

if strcmp(plots.density.type, 'contour')
    plots.density.range.on      = true;  % manual setting of ranges of the plots
    plots.density.range.x_min   = -3.5;
    plots.density.range.x_max   = 3.5;
    plots.density.range.y_min   = -3.5;
    plots.density.range.y_max   = 3.5;
end

plots.expect.energies.min   = -15;       % manual setting of minimum of energy expct. plot
plots.expect.energies.max   = +30;       % manual setting of maximum of energy expct. plot
