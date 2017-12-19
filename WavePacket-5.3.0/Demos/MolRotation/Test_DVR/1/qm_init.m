% Copyright (C) 2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (toggle)
global atomic hamilt plots psi space time

util.disp( '***********************************************' )
util.disp( 'Model system for testing spherical coordinate  ' )
util.disp( 'calculations. See model I in                   ' )
util.disp( ' J.Chem.Phys 95, 7392 (1991)                   ' )
util.disp( '***********************************************' )

% Only one angular degree of freedom. We set R to 1 and the mass such
% that hbar^2 / 2mR^2 = 1.47822 cm^-1
space.dof{1}       = grid.legendre;      % grid for expansion in Legendre polynomials
space.dof{1}.label = 'cos \Theta';
space.dof{1}.R_0   = 1;                  % constant value for R
space.dof{1}.m_0   = 0;                  % minor quantum number fixed to 0
space.dof{1}.l_max = 50;                 % maximum angular momentum
space.dof{1}.mass  = 1 / (2 * 1.47822 / atomic.w.cm_1); % adjusted mass

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 150;                   % Index of final time step

time.main.delta = 10000;                 % Size of big time steps: 240fs
time.sub.n      = 400;                   % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Operator splitting
time.propa.order = 3;                    % Strang method

% Hamiltonian operator 
hamilt.truncate.min = -2.0;              % Lower truncation of energy
hamilt.truncate.max = +1.5;              % Upper truncation of energy

hamilt.pot.handle    = @pot.metiu1;      % Potentials: Cylindrical coordinates
hamilt.pot.a0 =+700/atomic.w.cm_1;       % constant offset
hamilt.pot.a1 =+100/atomic.w.cm_1;       % first coefficient
hamilt.pot.a2 =-600/atomic.w.cm_1;       % second coefficient

% Initial wave function is a spherical harmonic of degree l=1.
psi.dof{1}.handle = @wav.fbr;            % use eigenstate of grid (spherical harmonic)
psi.dof{1}.state  = 2;                   % use second eigenstate (for m=0, this is l=1)

% Turn off plotting for speed.
plots.density.on = toggle;
plots.spectrum.on = false;
plots.expect.on = toggle;

% For qm_bound only
psi.eigen.stop = 15;
