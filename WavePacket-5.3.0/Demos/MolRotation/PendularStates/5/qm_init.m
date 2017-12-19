% Copyright (C) 2008-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init
global hamilt plots psi space

util.disp ( '***************************************' )
util.disp ( 'Pendular States for \Delta \omega = 80 ' )
util.disp ( '***************************************' )


% Note that the paper uses normalised units. We choose the rotational
% constant B = \hbar^2/(2mR^2) to be one a.u., so that the time units in
% the paper are the same numbers as the time units in a.u. used in the
% WavePacket run.

% Spatial discretization
space.dof{1} = grid.legendre;            % Legendre polynomials in cos theta
space.dof{1}.label = 'cos \Theta';
space.dof{1}.R_0 = 1;                    % constant value for R
space.dof{1}.m_0 = 0;                    % minor quantum number, fixed to 0
space.dof{1}.l_max = 201;                % maximum angular momentum/ number of points
space.dof{1}.mass = 0.5;                 % adjusted mass

% Orientation: cosine projector
space.amo{1}.handle = @amo.cosine;
space.amo{1}.exp = 1;

% Alignment: cosine^2 projector
space.amo{2}.handle = @amo.cosine;
space.amo{2}.exp = 2;

% Hamiltonian operator 
% hamilt.eigen.symmetry = 'u';             % may be used if potential symmetric
hamilt.pot.handle      = @pot.taylor;    % Taylor series in cos theta
hamilt.pot.v{1,1} = [0;-2*80];           % Force constant

% Select eigen/values/functions
psi.eigen.start        =  0;             % Lower index
psi.eigen.stop         = 15;             % Upper index

% Modify settings for appearance of plots (if desired)
plots.density.type = 'polar';            % polar plot

plots.density.pot.max = 200;             % customize density plot
plots.expect.energies.max = 300;         % Set maximum for energy plot
