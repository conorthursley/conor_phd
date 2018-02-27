% Copyright (C) 2008-2009 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time

util.disp ( '*******************************************************' )
util.disp ( 'Collinear model for the reaction F+D2 -> DF + D        ' )
util.disp ( 'see also Whitlock, Muckerman, JCP 61:4618 (1974) and   ' )
util.disp ( 'Muckerman, JCP 54:1155 (1971)                          ' )
util.disp ( '                                                       ' )
util.disp ( 'LEPS surface taken from Steinfeld, Francisco, Hase     ' )
util.disp ( '"Chemical Kinetics and Dynamics" p. 245, problem 7.14  ' )
util.disp ( '                                                       ' )
util.disp ( 'High collision energy                                  ' )
util.disp ( '*******************************************************' )

% atomic masses
m_F = atomic.mass.F19;                   % mass of fluorine (a.u.)
m_D = atomic.mass.H2;                    % mass of deuterium(a.u.)

% Spatial discretization
space.dof{1}       = grid.fft;           % D-D distance, equally spaced grid
space.dof{1}.label = 'D-D';
space.dof{1}.mass  = 1;                  % the internal kinetic energy is disabled
space.dof{1}.x_min = 0.4;                % minimum D-D distance
space.dof{1}.x_max = 7;                  % maximum D-D distance
space.dof{1}.n_pts = 128;                % number of points

space.dof{2}       = grid.fft;           % D-F distance, equally spaced grid
space.dof{2}.label = 'D-F';
space.dof{2}.mass  = 1;                  % internal KE is disabled anyway
space.dof{2}.x_min = 0.6;                % minimum D-F distance
space.dof{2}.x_max = 8;                  % maximum D-F distance
space.dof{2}.n_pts = 128;                % number of points

space.amo{1}.handle = @amo.reaction;     % population of reactant states
space.amo{1}.reac = 2;                   % educt distance is DOF #2
space.amo{1}.prod = 1;                   % product distance is DOf #1
space.amo{1}.side = 'r';                 % projection on the reactants
space.amo{1}.label= 'D_2 + F';           % labeling the reactants

space.amo{2}.handle = @amo.reaction;     % population of product states
space.amo{2}.reac = 2;                   % educt distance is DOF #2
space.amo{2}.prod = 1;                   % product distance is DOf #1
space.amo{2}.side = 'p';                 % projection on the products
space.amo{2}.label= 'DF + D';            % labeling the products

% Temporal discretization
time.main.start    = 0;                  % first time step
time.main.stop     = 100;                % last time step

time.main.delta    = 1/atomic.t.fs;      % length of one time step
time.sub.n         = 50;                 % number of sub steps

% Propagator
time.propa.handle = @ket.splitting;      % split operator 
time.propa.order = 3;                    % Strang splitting

% Hamiltonian operator 
hamilt.truncate.min= -0.3;               % lower truncation of energy
hamilt.truncate.max=  0.3;               % upper truncation of energy

hamilt.kin{1} = kin.triatomic;           % triatomic kinetic energy operator
hamilt.kin{1}.dof  = [1 2];              % act on the first two coordinates
hamilt.kin{1}.mass = [m_D m_D m_F];      % masses of the particles
hamilt.kin{1}.theta = pi;                % bending angle (pi for linear case)

hamilt.pot.handle  = @pot.leps;          % LEPS surface
hamilt.pot.theta = pi;                   % bending angle
hamilt.pot.d_e = [0.1745 0.2251 0.2251]; % diss. energy
hamilt.pot.r_e = [1.402  1.7329 1.7329]; % equilib. dist.
hamilt.pot.a   = [1.0277 1.1742 1.1742]; % range parameter
hamilt.pot.s   = [0.106  0.167  0.167];  % Sato parameter

% Absorbing boundary conditions
hamilt.nip.handle  = @nip.power;         % Absorbing boundary conditions
hamilt.nip.exp  = [2 2];
hamilt.nip.min  = [0.6 0.8];
hamilt.nip.max  = [6 7];

% Initial wave function
% We assume D-D in the HO ground state, and F coming in at a
% low speed (E_kin ~ 25 meV).
%
% Note that the coordinate transformation also affects the
% initial wave function!
%
% In atomic positions (x_a = Fluorine, x_b = first deuterium, x_c = second deuterium)
% the initial wave function is simply
% 
% exp(i k_0 x_a) * exp(-i k_0 0.5*(x_b+x_c) * exp(-a(1/2(x_b-x_c) - R0)^2)
%                * f(x_a) * g(0.5(x_b + x_c))
%
% Read: D2 and F crash into each other with momentum k_0 each (CMS frame), 
% D2 is in the vibrational ground state, and F and D2 are located somewhere
% (described by the functions f and g). If we transform this to bond length
% coordinates, we find that the first two terms give
%
% exp(-i k_0 R_ab) * exp(-i k_0/2 R_bc),
%
% i.e., the initial wave function for the D-D distance needs to get half
% the momentum, too. The transformations of f and g are nontrivial, too.
% Here, we just assume a Gaussian for the D-F distance.

omega = sqrt(2*0.174*1.0277^2/(m_D/2));

psi.dof{1}.handle = @wav.gauss;
psi.dof{1}.pos_0  = 1.4019;
psi.dof{1}.mom_0  = -3.5;
psi.dof{1}.width  = sqrt(1/(2*(m_D/2)*omega));

psi.dof{2}.handle = @wav.gauss;
psi.dof{2}.pos_0  = 5;
psi.dof{2}.mom_0  = -7;
psi.dof{2}.width  = 0.3;

% Modify settings for appearance of plots (if desired)
plots.density.type = 'surface'; 
plots.density.surface.view = [110 70];
plots.expect.energies.max = 0.2;
