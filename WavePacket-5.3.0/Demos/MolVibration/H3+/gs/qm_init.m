% Copyright (C) 2008 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time

util.disp( '**********************************************' )
util.disp( 'Testing case for triatomic Jacobi Hamiltonian ' )
util.disp( 'Lowest state of H3+                           ' )
util.disp( '**********************************************' )

% First coordinate: H-H distance
space.dof{1}       = grid.fft;           % equally spaced FFT grid
space.dof{1}.label = 'H-H';              % label
space.dof{1}.mass  =atomic.mass.H1/2;    % reduced mass of coordinate
space.dof{1}.x_min = 0.9;                % minimum grid point
space.dof{1}.x_max = 4.5;                % maximum grid point
space.dof{1}.n_pts = 64;                 % number of grid points

% Second coordinate: Distance between H and H2 center (diss. coordinate)
space.dof{2}       = grid.fft;           % equally spaced FFT grid
space.dof{2}.label = 'H-H2';             % label
space.dof{2}.mass  = 2/3*atomic.mass.H1; % reduced mass of coordinate
space.dof{2}.x_min = 0.7;                % minimum grid point
space.dof{2}.x_max = 4;                  % maximum grid point
space.dof{2}.n_pts = 64;                 % number of grid points

% Third coordinate: angle (/cos of angle) between the former two
space.dof{3}       = grid.legendre;      % Legendre/Spherical harmonics expansion of cos theta
space.dof{3}.label = 'cos \Theta';       % label
space.dof{3}.mass  = 1;                  % unused
space.dof{3}.m_0   = 0;                  % minor quantum number of spherical harmonic
space.dof{3}.R_0   = 1;                  % unused
space.dof{3}.l_max = 40;                 % maximum quantum number

% Temporal discretization
time.main.start = 00;                    % Index of initial time step
time.main.stop  = 10;                    % Index of final time step

time.main.delta = 5/atomic.t.fs;         % Size of main time step: 5fs
time.sub.n      = 1;                     % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.cheby_imag;     % Chebyshev expansion in imaginary time
time.propa.order = 0;                    % auto-determine number of polynomials
time.propa.precision = 1e-10;            % requested precision of the propagation

% Hamiltonian operator 
hamilt.truncate.min = 0;                 % Lower truncation of energy
hamilt.truncate.max = +1;                % Upper truncation of energy

hamilt.pot.handle = @pot.H3plus;         % Use H3+ potential
hamilt.pot.c_dof = 3;                    % 3rd DOF is angular coordinate

hamilt.kin{1} = kin.jacobi;              % Use kin.jacobi for the angular coordinate
hamilt.kin{1}.dof_r = 1;                 % first DOF is bonding distance
hamilt.kin{1}.dof_R = 2;                 % second DOF is dissociative one
hamilt.kin{1}.dof_c = 3;                 % third DOF is cos(angle)
hamilt.kin{1}.mass_r = atomic.mass.H1/2; % reduced mass of the bonding DOF
hamilt.kin{1}.mass_R =2/3*atomic.mass.H1;% reduced mass of the dissociative DOF

% Initial wave function is a Gaussian. It is located at the
% equilibrium positions of the potential. The angular coordinate
% is shifted so that we avoid a symmetric setup
psi.corr.handle = @wav.gauss;
psi.corr.pos_0  = [1.6 1.4 cosd(80)];
psi.corr.mom_0  = [0 0 0];
psi.corr.width  = [0.6 0.6 0.2];
psi.corr.factor = [1 1 1];

% We do not have any reasonable plotting facilities, so turn it off.
plots.density.on = false;
plots.expect.on  = false;
