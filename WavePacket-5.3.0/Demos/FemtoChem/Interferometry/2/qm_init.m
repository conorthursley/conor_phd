% Copyright (C) 2009,2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global atomic hamilt plots psi space time

util.disp ('****************************************')
util.disp ('Calculation of the pump only signal for ')
util.disp ('excitation of Na2.                      ')
util.disp ('see J.Chem.Phys. 100:5448               ')
util.disp ('****************************************')

% Number of coupled equations
hamilt.coupling.n_eqs  = 3;
hamilt.coupling.representation = 'dia';
hamilt.coupling.labels = {'X^1\Sigma_g^+', 'A^1\Sigma_u^+', '2^1\Pi_g'};

% Grid definition
space.dof{1}       = grid.fft;           % using FFT grid
space.dof{1}.n_pts = 256;                % number of points
space.dof{1}.x_min = 4;                  % lower boundary of the grid
space.dof{1}.x_max = 14;                 % upper boundary of the grid
space.dof{1}.mass  = atomic.mass.Na23/2; % reduced mass of 23Na2

% Temporal discretisation
time.main.start = 0;                     % index of first time step
time.main.stop  = 210;                   % index of last time step

time.main.delta = 10/atomic.t.fs;        % 10 fs per time step
time.sub.n      = 500;                   % propagation time step 20 as

% Propagator
time.propa.handle = @ket.splitting;      % Split operator
time.propa.order = 3;                    % Strang splitting

% Electric field is a sin^2 pulse 
time.efield.shape   = 'gauss';           % a Gauss pulse
time.efield.fwhm    = 30/atomic.t.fs;    % with 30 fs FWHM
time.efield.delay   = 100/atomic.t.fs;   % delayed by 100 fs
time.efield.polar   = 0;                 % fixed orientation => irrelevant
time.efield.phase   = 0;                 % do not bother with phase shift
time.efield.frequ   = 0.073;             % frequency (corr. to 625 nm)
time.efield.ampli   = sqrt(10/atomic.I.GW_cm2);% 10 GW/cm^2

% Hamiltonian operator
hamilt.truncate.min = -0.03;             % lower truncation of energy
hamilt.truncate.max =  0.25;             % upper truncation of energy

hamilt.pot.handle = @pot.interp;         % potential energy
hamilt.pot.pos_conv = atomic.l.A;        % conversion factor for coordinates
hamilt.pot.pot_conv = atomic.E.eV;       % conversion factor for energies

hamilt.dip.handle = @FC_dip;             % dipole moments; see below.

% Initial wave function; Interpolate from input file.
psi.corr.handle = @wav.interp; 

% Modify appearance of plots
plots.density.type        = 'contour';   % Draw contour lines

plots.density.contour.min = -0.01;       % bounds for the contour lines
plots.density.contour.max = 0.01;

plots.expect.energies.max = 0.17;


% Save the wave function
psi.save.export = true;
psi.save.dir    = '.';
psi.save.file   = 'pump';
% if the output is really huge (long time series or so), you might
% want to break it into smaller chunks, which also speeds up the loading
% afterwards (because smaller files are loaded into memory)
%psi.save.mem = 10 * 2^20;      % new file every 10 MB data

end                                      % qm_init


% Now specify a function that gives a simple Franck-Condon transition for
% the different states. Unfortunately, the reference does not give values here,
% so we just invent some vaguely based on reality.
function FC_dip
    global hamilt space
    hamilt.d_x.grid_ND{1,1} = [];
    hamilt.d_x.grid_ND{2,2} = [];
    hamilt.d_x.grid_ND{3,3} = [];
    hamilt.d_x.grid_ND{1,2} = ones(size(space.dvr.grid_ND{1}));
    hamilt.d_x.grid_ND{2,3} = ones(size(space.dvr.grid_ND{1}));
end

