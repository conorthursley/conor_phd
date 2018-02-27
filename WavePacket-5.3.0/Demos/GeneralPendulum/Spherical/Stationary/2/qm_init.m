% Copyright (C) 2016 Burkhard Schmidt
%
% This file is under public domain. You may modify or incorporate it into other
% works without any restrictions.

function qm_init (b,m)
global hamilt plots space

util.disp ( '*********************************************' )
util.disp ( 'Supersymmetry and eigensurface topology' )
util.disp ( 'of the spherical quantum pendulum' )
util.disp ( 'B. Schmidt and B. Friedrich' )
util.disp ( 'Phys. Rev. A 91, 022111' )
util.disp ( 'DOI:10.1103/PhysRev.A.91.022111' )
util.disp ( 'Reproducing red circles in Fig. 3' )
util.disp ( 'USING HERE: Gauss-Legendre-DVR for COS(THETA)' )
util.disp ( '*********************************************' )

% Spatial discretization
space.dof{1} = grid.legendre;            % Gauss-Legendre DVR in cos(theta)
space.dof{1}.label = 'cos \Theta';
space.dof{1}.R_0 = 1;                    % constant value for R
space.dof{1}.m_0 = m;                    % minor quantum number 
space.dof{1}.l_max = 100;                % maximum angular momentum/ number of points
space.dof{1}.mass = 0.5;                 % adjusted mass

% Orientation: cosine projector
space.amo{1}.handle = @amo.cosine;
space.amo{1}.exp = 1;

% Alignment: cosine^2 projector
space.amo{2}.handle = @amo.cosine;
space.amo{2}.exp = 2;

% Pendular potential
hamilt.pot.handle = @pot.taylor;         % Taylor series in cos(theta)
hamilt.pot.v {1,1} = [-2*b*(m+1);-2*b^2];% eta and zeta parameters

% Modify settings for appearance of plots (if desired)
plots.density.type = 'curve';            % polar plot
plots.density.pot.max = 200;             % customize density plot
plots.expect.population.min = -0.3;      % customize population/amo plot
plots.expect.population.max = +1.1;      % customize population/amo plot
plots.expect.energies.min = -150;        % customize energy plot
plots.expect.energies.max = +150;        % customize energy plot
