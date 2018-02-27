% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '***********************' )
util.disp ( 'Schroedingers cat state' )
util.disp ( '***********************' )

% Spatial discretization
space.dof{1} = grid.fft;                 % using fft grid
space.dof{1}.mass = 1;                   % mass
space.dof{1}.n_pts = 128;                % Number of grid points
space.dof{1}.x_min = -15.0;              % Lower bound of grid 
space.dof{1}.x_max = +15.0;              % Upper bound of grid

% Temporal discretization band propagator
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 010;                   % Index of final time step

time.main.delta  = pi/10;                % Size of time steps 
time.sub.n   =  100;                     % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    =   0.0;          % Lower truncation of energy
hamilt.truncate.max    =  30.0;          % Upper truncation of energy

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width = [  1; 1 ];            % Width 
psi.dof{1}.pos_0 = [  0; 0 ];            % Center in position representation
psi.dof{1}.mom_0 = [ -2; 2 ];            % Center in momentum representation

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % Contour plot

plots.density.range.on = true;           % Manual setting of plotting range
plots.density.range.x_min = -12;
plots.density.range.x_max = 12;
plots.density.range.y_min = -4;
plots.density.range.y_max = 4;

plots.expect.energies.max = 2.5;
