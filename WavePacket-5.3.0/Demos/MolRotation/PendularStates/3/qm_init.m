% Copyright (C) 2008 Ulf Lorenz, Burkhard Schmidt
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init(Delta)
global hamilt plots psi space time

util.disp( '***********************************************' )
util.disp( 'Time evolution of pendular states              ' )
util.disp( '\sigma = 5 pulse and varying intensity         ' )
util.disp( ' J.Chem.Phys 110, 3870 (1999)                  ' )
util.disp( '***********************************************' )


% Note that the paper uses normalised units. We choose the rotational
% constant B = \hbar^2/(2mR^2) to be one a.u., so that the time units in
% the paper are the same numbers as the time units in a.u. used in the
% WavePacket run.

% One angular degree of freedom
space.dof{1} = grid.legendre;            % Gauss-Legendre DVR in cos theta
space.dof{1}.label = 'cos \Theta';
space.dof{1}.R_0 = 1;                    % constant value for R
space.dof{1}.m_0 = 0;                    % minor quantum number fixed to 0
space.dof{1}.l_max = 10;                 % maximum angular momentum
space.dof{1}.mass = 0.5;                 % adjusted mass

% Alignment: cosine^2 projector
space.amo{1}.handle = @amo.cosine;
space.amo{1}.exp = 2;

% Temporal discretization
time.main.start = 000;                   % Index of initial time step
time.main.stop  = 500;                   % Index of final time step

time.main.delta = 0.1;                   % Size of big time steps
time.sub.n      = 100;                   % Number of sub steps per time step

% Propagator
time.propa.handle = @ket.splitting;      % split-operator method
time.propa.order = 3;

% Electric fild with properly designed parameters
time.efield.polar = 0;                   % Required by the propagator
time.efield.shape = 'gauss';             % Gaussian pulse
time.efield.delay = 15;                  % pulse delay, measured from paper
time.efield.fwhm  = 5*sqrt(8*log(2));    % FWHM of pulse
time.efield.frequ = 0;                   % no oscillations
time.efield.phase = 0;                   % no phase shift
time.efield.ampli = sqrt(Delta);         % amplitude of field

% Hamiltonian operator: polarizability interaction 
hamilt.pol.handle = @pol.taylor;         % Taylor series in cos Theta
hamilt.pol.x.v{1,1} = [0;2];             % Set the polarizability to cos^2

% Initial wave function is a spherical harmonic of degree l=1.
psi.dof{1}.handle = @wav.fbr;            % use eigenstate of grid (spherical harmonic)
psi.dof{1}.state  = 1;                   % use first eigenstate (for m=0, this is l=0)

% Plot nothing. This is done in a summary at the end of the runs.
plots.density.on  = false;
plots.spectrum.on = false;
plots.expect.on   = false;
