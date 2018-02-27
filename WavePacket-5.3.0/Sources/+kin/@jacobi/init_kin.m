% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008,2011 Ulf Lorenz
%
% see the README file for license details.

function obj = init_kin (obj, fraction, output)

% initialises the kinetic energy

global hamilt space time

if nargin < 3
    output = true;
end

%% Checks

if isempty( obj.mass_R ) || isempty( obj.mass_r )
     util.error ('masses not set')
end
if isempty( obj.dof_c ) 
     util.error ('c degree of freedom not set')
end
if isempty( obj.dof_R ) 
     util.error ('R degree of freedom not set')
end
if isempty( obj.dof_r ) && isempty( obj.r_0 )
     util.error ('r degree of freedom not set')
end

% Check that we have a reasonable grid (i.e. Legendre polynomials)
if ~isa ( space.dof{obj.dof_c}, 'grid.legendre' )
     util.error ('Jacobi kinetic energy only works for Legendre grids')
end


%% Informational output

if output
    util.disp (' ')
    util.disp ('*******************************************************')
    util.disp ('Kinetic energy: Bending c of triatomic molecule ABC    ')
    util.disp ('                                                       ')
    util.disp ('           [     1             1    ] ^ 2              ')
    util.disp (' T (c) = - [ --------  +   -------- ] L                ')
    util.disp ('           [ 2 M R^2       2 m r^2  ]                  ')
    util.disp ('                                                       ')
    util.disp ('where R, M is the distance and reduced mass of B and C,')
    util.disp ('r distance of A to CMS of B and C, m the reduced mass  ')
    util.disp ('of A and BC. See for example J. Chem. Phys 116:4403    ')
    util.disp ('*******************************************************')
    util.disp ( [ 'M     : ' num2str(obj.mass_R) ] )
    util.disp ( [ 'm     : ' num2str(obj.mass_r) ] )
    util.disp ( [ 'DOF c : ' num2str(obj.dof_c)  ] )
    util.disp ( [ 'DOF R : ' num2str(obj.dof_R)  ] )
    if ~isempty(obj.dof_r)
        util.disp ( [ 'DOF r : ' num2str(obj.dof_r) ] )
    else
        util.disp ( [ 'constant r : ' num2str(obj.r_0) ] )
    end
end


%% Create all the grids

% Prefactor
if isempty( obj.r_0 )
    obj.grid = 1 ./ (2 * space.dvr.grid_ND{obj.dof_R}.^2 * obj.mass_R) ...
               + 1 ./ (2 * space.dvr.grid_ND{obj.dof_r}.^2 * obj.mass_r);
else
    obj.grid = 1 ./ (2 * space.dvr.grid_ND{obj.dof_R}.^2 * obj.mass_R) ...
               + 1 ./ (2 * obj.r_0.^2 * obj.mass_r);
end

% L^2
obj.grid = space.fbr.grid_ND{obj.dof_c}.^2 .* obj.grid;

% Truncation
obj.grid(obj.grid > hamilt.truncate.delta) = hamilt.truncate.delta;

% Short-time propagator
obj.grid_exp = exp(-1i * obj.grid * time.sub.delta * fraction);
