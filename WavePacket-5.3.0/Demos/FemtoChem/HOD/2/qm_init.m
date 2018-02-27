% Copyright (C) 2010 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time params


% This file only works with the associated 'rundemo.m' script.  It uses a
% global variable "params" from which the time stepping, wave function
% initialisation, and wave function saving behaviour are copied. The rundemo.m
% script then recycles this initialise.m function to send the initial and target
% wavepackets back and forth in time.

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

% Temporal discretisation copied from the params global variable.
% The electric field is also copied.
time = params.time;

% Propagator
time.propa.handle = @ket.splitting;      % split operator method
time.propa.order = 3;                    % use Strang splitting

% Hamiltonian operator
hamilt.truncate.min = -0.5;              % lower truncation of energy
hamilt.truncate.max = 0.5;               % upper truncation of energy

hamilt.kin{1} = kin.triatomic;           % Kinetic energy for fixed bending angle
hamilt.kin{1}.theta = 104.52 * pi/180;   % bending angle
hamilt.kin{1}.dof = [1 2];               % indices of coordinates for AB, BC distance
hamilt.kin{1}.mass= [2.13 16.00 1.00]/ atomic.m.u;
                                         % note that the D mass is off by a few percent
                                         % due to approximations in Ashwanis setup
                                         % (he assumed mu_OD = 2 * mu_OH).

hamilt.pot.handle = @pot.H2O;            % potential energy

hamilt.dip.handle = @dip.H2O;            % transition dipole moments

% Initial wave functionand some details on whether or not we save the wave
% function are also copied from the externally setup params variable.
psi = params.psi;

% Typically, we do not want to have plots, but let the rundemo script decide.
% If it does not prohibit plotting, there are a few standard settings we still
% can do.
if params.plot == false
	plots.density.on = false;
	plots.expect.on = false;
end

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
