% Copyright (C) 2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (toggle)
global atomic hamilt plots psi space time

util.disp( '***********************************************' )
util.disp( 'Model system for testing spherical coordinate  ' )
util.disp( 'calculations. See model II of HCl in           ' )
util.disp( ' J.Chem.Phys 95, 7392 (1991)                   ' )
util.disp( '***********************************************' )

% Reduced mass
mass_HCl = 0.97958 / atomic.m.u;         

% Grid parameters for HCl problem
space.dof{1}       = grid.fft;           % Radial grid
space.dof{1}.label = 'H-Cl';             % Label
space.dof{1}.mass  = mass_HCl;           % Mass of kinetic energy operator
space.dof{1}.x_min = 0.85 / atomic.l.A;  % grid start
space.dof{1}.x_max = 4.5  / atomic.l.A;  % grid end
space.dof{1}.n_pts = 96;                 % number of grid points

space.dof{2}       = grid.legendre;      % angular grid
space.dof{2}.label = 'cos \Theta';
space.dof{2}.R_dof = 1;                  % couple to first DOF in kinetic energy
space.dof{2}.m_0   = 0;                  % minor quantum number fixed to 0
space.dof{2}.l_max = 30;                 % maximum angular momentum
space.dof{2}.mass  = mass_HCl;           % mass of kinetic energy operator

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 120;                   % Index of final time step

time.main.delta = 5000;                  % Size of big time steps: 120 fs
time.sub.n      = 1000;                  % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % Operator splitting
time.propa.order = 3;                    % Strang method

% Hamiltonian operator 
hamilt.truncate.min = 0;                 % Lower truncation of energy
hamilt.truncate.max = 1;                 % Upper truncation of energy

hamilt.pot.handle     = @pot.morse;      % Morse potential
hamilt.pot.d_e = [37244/atomic.w.cm_1; 0];   % dissociation energy
hamilt.pot.r_e = [1.2746/atomic.l.A; 0]; % equilibrium length
hamilt.pot.alf = [2.38/(1.2746/atomic.l.A); 0];  % range parameter
hamilt.pot.t_e = 0;               % energy shift

% Initial wave function is a Gaussian in R and a superposition of
% Spherical harmonics in Theta
psi.dof{1}.handle = @wav.gauss;          % Gaussian
psi.dof{1}.pos_0  = 1.7746/atomic.l.A;   % mean position
psi.dof{1}.mom_0  = 0;                   % mean momentum
psi.dof{1}.width  = 1/sqrt(2*mass_HCl ...
                            *2989.78/atomic.w.cm_1);% width of Gaussian

psi.dof{2}.handle = @wav.fbr;            % Superposition of eigenstates
psi.dof{2}.coeffs = [ 0         0.261794  0.159839 ...
                          -0.181825 -0.295285 -0.043585 ...
                           0.260620  0.271072 -0.005645 ...
                          -0.277002 -0.302614 -0.090726 ...
                           0.182073  0.355587  0.383720 ...
                           0.311739  0.208091  0.118826 ...
                           0.059407  0.026403  0.010545 ...
                           0.003817  0.001260  0.000381 ...
                           0.000106  0.000027  0.000007];

% Turn off plotting for speed.
plots.density.on = toggle;
plots.spectrum.on = false;
plots.expect.on = toggle;

% For qm_bound only
psi.eigen.stop = 200;

