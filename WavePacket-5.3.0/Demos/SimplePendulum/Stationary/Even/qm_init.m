% Copyright (C) 2004-2007 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '***************************************' )
util.disp ( 'Eigenstates of plane quantum pendulum  ' )
util.disp ( '***************************************' )

% Spatial discretization
space.dof{1}       = grid.fft;           % using fft grid
space.dof{1}.mass  = 1;                  % mass
space.dof{1}.n_pts = 192;                % Number of grid points
space.dof{1}.x_min = 0;                  % Lower bound of grid 
space.dof{1}.x_max = 2*pi;               % Upper bound of grid

space.amo{1}.handle = @amo.cosine;       % Degree of orientation

% Hamiltonian operator 
hamilt.truncate.min  = 000;              % Lower truncation of energy
hamilt.truncate.max  = 300;              % Upper truncation of energy

hamilt.pot.handle    = @pot.pendulum;    % Intramolecular torsion
hamilt.pot.eta = -50;                    % Prefactor of cos-potential
hamilt.pot.v_0 = +50;                    % Energy offset

% Fourier grid Hamiltonian
hamilt.eigen.symmetry  = 'g';            % Even parity states

% Select eigen/values/functions
psi.eigen.start        =  0;             % Lower index
psi.eigen.stop         = 20;             % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type        = 'polar';
plots.expect.energies.max = 300;
plots.density.export.on = true;
plots.expect.export.on = true;
