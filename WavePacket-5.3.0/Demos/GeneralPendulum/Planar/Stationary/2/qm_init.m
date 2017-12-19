% Copyright (C) 2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '***************************************' )
util.disp ( 'Supersymmetry and eigensurface topology' )
util.disp ( 'of the planar quantum pendulum' )
util.disp ( 'B. Schmidt and B. Friedrich' )
util.disp ( 'Front. Phys. 2, 37' )
util.disp ( 'DOI:10.3389/fphy.2014.00037' )
util.disp ( 'Reproducing red circles in Fig. 4' )
util.disp ( 'USING HERE: FFT-DVR for THETA' )
util.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % Fourier grid
space.dof{1}.n_pts = 256;                % Number of grid points
space.dof{1}.x_min = -2*pi;              % Lower bound of grid
space.dof{1}.x_max = +2*pi;              % Upper bound of grid
space.dof{1}.mass  =  1/2;               % Mass for the kinetic energy

% Orientation
space.amo{1}.handle = @amo.cosine;       % cosine projector
space.amo{1}.exp = 1;                    % exponent

% Alignment
space.amo{2}.handle = @amo.cosine;       % cosine^2 projector
space.amo{2}.exp = 2;                    % exponent

% Potential: beta=5, kappa=6
hamilt.pot.handle = @pot.pendulum;       % potential for generalized pendula
hamilt.pot.eta  = 30;                    % Orientation: cos
hamilt.pot.zeta = 25;                    % Alignment: cos^2

% Generate first 25 eigenstates
psi.eigen.start = 0;
psi.eigen.stop = 25;

% Modify settings for appearance of plots (if desired)
plots.density.type = 'contour';          % Wigner contour plot
plots.density.pot.max = 050;             % customize density plot
plots.expect.energies.max = 040;         % Set maximum for energy plot
