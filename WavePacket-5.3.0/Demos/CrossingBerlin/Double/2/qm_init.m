% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space time

util.disp ( '**********************************************' )
util.disp ( 'Dual crossing example (Horenko, Schmidt: 2004)' )
util.disp ( '   ' )
util.disp ( '             ( A     C   ) ' )
util.disp ( ' V_dia (R) = (           ) ' )
util.disp ( '             ( C    BR^2 ) ' )
util.disp ( '   ' )
util.disp ( 'with A=0.1, B=0.05, C=0.02' )
util.disp ( '   ' )
util.disp ( '**********************************************' )

% Number of (coupled) Schrödinger equations and representation
hamilt.coupling.n_eqs = 2;
hamilt.coupling.representation = 'adi';

% Spatial discretization
space.dof{1} = grid.fft;                 % using FFT grid
space.dof{1}.mass = 2000;                % mass
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -5.0;               % Lower bound of grid 
space.dof{1}.x_max =  5.0;               % Upper bound of grid

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 060;                   % Index of final time step

time.main.delta  = 10;                   % Size of time steps 
time.sub.n       = 050;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Split operator method
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min    = 0.00;           % Lower truncation of energy
hamilt.truncate.max    = 1.00;           % Upper truncation of energy

hamilt.pot.handle = @pot.taylor;         % Taylor series: for dual crossing
hamilt.pot.c = cell(2);                  % Pre-allocate cell array
hamilt.pot.v = cell(2);                  % Pre-allocate cell array
hamilt.pot.c{1,1} = 0.1;                 % Constant potential
hamilt.pot.v{2,2} = [0;2*0.05];          % Parabolic potential
hamilt.pot.c{1,2} = 0.02;                % Coupling (constant)

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % Negative imaginary potential
hamilt.nip.exp  = 4;                     % Exponent
hamilt.nip.min = -4.0;                   % Beginning of inner grid region
hamilt.nip.max = +4.0;                   % End of inner grid region

% Initial wave function
psi.dof{1}.handle= @wav.gauss;           % Gaussian-shaped wavepacket
psi.dof{1}.width =  0.2991/2;            % Width 
psi.dof{1}.pos_0 = -3.0;                 % Center in position representation
psi.dof{1}.mom_0 = 20.0;                 % Center in momentum representation

psi.init.representation = 'adi';
psi.init.coeffs       = [1 0];           % Initially only lower adiabatic state populated

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';

plots.density.surface.view  = [-5 25];   % View point for surface plot: [az el]

if strcmp(plots.density.type, 'contour')
    plots.density.range.on    = true;    % manual setting of plotting range
    plots.density.range.x_min = -4.5;
    plots.density.range.x_max = 4.5;
    plots.density.range.y_min = -30;
    plots.density.range.y_max = 35;
end

plots.expect.energies.max = 0.4;
