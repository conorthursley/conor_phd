% Copyright (C) 2009 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time


util.disp ('********************************')
util.disp ('Bond breakage of HOD            ')
util.disp ('with an unshaped laser pulse.   ')
util.disp ('see PhysRev A 78:065402 (2008)  ')
util.disp ('********************************')

% Number of coupled equations
hamilt.coupling.n_eqs  = 2;
hamilt.coupling.representation = 'dia';
hamilt.coupling.labels = {'X^1 A_1', 'A^1 B_1'};

% Grid definition; masses are not needed here since we use
% a custom operator; it would be smarter to use an unbalanced
% grid (O-D is heavier and thus needs more grid points than O-H),
% but we shall stick to the original setup as close as possible.
space.dof{1}       = grid.fft;           % using FFT grid
space.dof{1}.label = 'O-D';              % label for axis descriptions
space.dof{1}.n_pts = 512;                % number of points
space.dof{1}.x_min = 1.15;               % lower boundary of the grid
space.dof{1}.x_max = 26.75;              % upper boundary of the grid

space.dof{2}       = grid.fft;           % using FFT grid
space.dof{2}.label = 'O-H';              % label for axis descriptions
space.dof{2}.n_pts = 512;                % number of points
space.dof{2}.x_min = 1.15;               % lower boundary of the grid
space.dof{2}.x_max = 26.75;              % upper boundary of the grid

space.amo{1}.handle = @amo.reaction;     % projecting on D + OH channel
space.amo{1}.reac  = 1;
space.amo{1}.prod  = 2;
space.amo{1}.side  = 'r';
space.amo{1}.label = 'D + OH';

space.amo{2}.handle = @amo.reaction;     % projecting on H + OD channel
space.amo{2}.reac  = 1;
space.amo{2}.prod  = 2;
space.amo{2}.side  = 'p';
space.amo{2}.label = 'H + OD';

% Temporal discretisation
time.main.start = 0;                     % index of first time step
time.main.stop  = 100;                   % index of last time step

time.main.delta = 20;                    % approx 0.5 fs per time step
time.sub.n      = 20;                    % propagation time step 1 a.u.

% Propagator
time.propa.handle = @ket.splitting;      % split operator method
time.propa.order = 3;                    % use Strang splitting

% Electric field is a sin^2 pulse 
time.efield.shape   = 'sin^2';           % sin^2 shape
time.efield.fwhm    = 25/atomic.t.fs;    % 25 fs  FWHM
time.efield.delay   = 25/atomic.t.fs;    % centered around t = 25 fs
time.efield.polar   = 0;                 % polarization of pulse
time.efield.ampli   = sqrt(18/atomic.I.TW_cm2);   % 18 TW/cm^2
time.efield.frequ   = 0.27;              % frequency
time.efield.phase   = 0;                 % phase shift

% Hamiltonian operator
hamilt.truncate.min = -0.5;              % lower truncation of energy
hamilt.truncate.max = +0.5;              % upper truncation of energy

hamilt.kin{1} = kin.triatomic;           % Kinetic energy for fixed bending angle
hamilt.kin{1}.theta = 104.52 * pi/180;   % bending angle
hamilt.kin{1}.dof = [1 2];               % indices of coordinates for AB, BC distance
hamilt.kin{1}.mass= [2.13 16.00 1.00]/ atomic.m.u;
                                         % note that the D mass is off by a few percent
                                         % due to approximations in Ashwanis setup
                                         % (he assumed mu_OD = 2 * mu_OH).

hamilt.pot.handle = @pot.H2O;            % potential energy

hamilt.dip.handle = @dip.H2O;            % transition dipole moments

% Absorbing boundary conditions
hamilt.nip.handle = @nip.power;          % negative imaginary potential
hamilt.nip.exp = [2 2];                  % Exponent
hamilt.nip.min = [1.3 1.3];              % Begin of inner grid region
hamilt.nip.max = [25 25];                % End of inner grid region

% Initial wave function; HOD ground state modelled by two Gaussians
psi.dof{1}.handle = @wav.gauss;          % Gaussian initial wave function
psi.dof{1}.width  = 1/sqrt(4*23.18);     % width of Gaussian
psi.dof{1}.pos_0  = 1.806;               % center in position representation
psi.dof{1}.mom_0  = 0;                   % center in momentum representation

psi.dof{2}.handle = @wav.gauss;          % Gaussian initial wave function
psi.dof{2}.width  = 1/sqrt(4*18.84);     % Width of Gaussian
psi.dof{2}.pos_0  = 1.833;               % Center in position representation
psi.dof{2}.mom_0  = 0;                   % center in momentum representation

psi.init.coeffs        = [1 0];          % Population of initial state
psi.init.representation= 'dia';          % Representation of initial wave function

% Modify appearance of plots
plots.density.type        = 'contour';   % Draw contour lines

plots.density.range.on    = true;        % manually set ranges
plots.density.range.x_min = 1.2;
plots.density.range.x_max = 6;
plots.density.range.y_min = 1.2;
plots.density.range.y_max = 6;

plots.density.contour.nlev = [50 15];    % number of contour lines
plots.density.rho_max.dvr = 6;           % factor that determines ranges of density
                                         % for contour lines to draw
plots.density.pot.min = -0.35;           % finetuning for drawing of energy curves
plots.density.pot.max = 0;

plots.expect.energies.max = 0.3;         % manually set ranges of energy plot
