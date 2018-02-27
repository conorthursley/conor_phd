% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '**********************************************' )
util.disp ( 'Razavy symmetric double well potential' )
util.disp ( 'Reference calculation using Gauss-Hermite DVR ' )
util.disp ( '**********************************************' )

% Spatial discretization
space.dof{1}       = grid.hermite;       % Use Gauss-Hermite grid
space.dof{1}.mass  = 1/2;                % Particle mass
space.dof{1}.n_pts = 64;                 % Number of grid points
space.dof{1}.v_2   = 7;                  % Frequency of oscillator

% Hamiltonian operator 
% hamilt.eigen.symmetry = 'g';             % Symmetry may be exploited
hamilt.truncate.min    = -15.0;          % Lower truncation of energy
hamilt.truncate.max    = 1000.0;         % Upper truncation of energy

% Razavy Potential: beta=0.1, kappa=-7
hamilt.pot.handle      = @pot.razavy;    % Hyperbolic potential
hamilt.pot.modified = true;
hamilt.pot.eta = -0.7;
hamilt.pot.zeta = 0.01;

% Select eigen/values/functions
psi.eigen.start        = 0;              % Lower index
psi.eigen.stop         = 30;             % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type = 'curve';            % simple plot of wavefunction
plots.density.pot.max = 200;
plots.expect.energies.max = 80;
