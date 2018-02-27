% Copyright (C) 2016 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (b,m)
global hamilt plots space

util.disp ( '***************************************' )
util.disp ( 'Supersymmetry and eigensurface topology' )
util.disp ( 'of the spherical quantum pendulum' )
util.disp ( 'B. Schmidt and B. Friedrich' )
util.disp ( 'Phys. Rev. A 91, 022111' )
util.disp ( 'DOI:10.1103/PhysRev.A.91.022111' )
util.disp ( 'Reproducing red circles in Fig. 3' )
util.disp ( 'USING HERE: FFT-DVR for THETA' )
util.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % Fourier grid
space.dof{1}.label = '\Theta';
space.dof{1}.n_pts = 1024;               % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid
space.dof{1}.x_max = +1*pi;              % Upper bound of grid
space.dof{1}.mass  =  1/2;               % Mass for the kinetic energy

% Orientation: cosine projector
space.amo{1}.handle = @amo.cosine;
space.amo{1}.exp = 1;

% Alignment: cosine^2 projector
space.amo{2}.handle = @amo.cosine;
space.amo{2}.exp = 2;

% Pendular potential
hamilt.pot.handle = @pot.pendulum;       % potential for generalized pendula
hamilt.pot.xi = m^2 - 1/4;               % Azimuthal rotation
hamilt.pot.eta = 2*b*(m+1);              % Orientation: cos Theta
hamilt.pot.zeta = b^2;                   % Alignment: cos^2 Theta
hamilt.pot.v_0 = -1/4;                   % Energy shift

% Modify settings for appearance of plots (if desired)
plots.density.type = 'curve';            % polar plot
plots.density.pot.max = 200;             % customize density plot
plots.expect.population.min = -0.3;      % customize population/amo plot
plots.expect.population.max = +1.1;      % customize population/amo plot
plots.expect.energies.min = -150;        % customize energy plot
plots.expect.energies.max = +150;        % customize energy plot

