% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2011 Ulf Lorenz
%               2012 Ulf Lorenz
%
% see the README file for license details.

function obj = init_kin (obj, fraction, output)

% initialises the kinetic energy

global hamilt space time


if nargin < 3
    output = true;
end


% output some information
if output
    util.disp (' ')
    util.disp ('*******************************************************')
    util.disp ('Kinetic energy: Triatomic ABC molecule (fixed angle)   ')
    util.disp ('                                                       ')
    util.disp ('                1    (   d    ) 2     1    (   d    ) 2')
    util.disp (' T (d/dR) = - ------ ( ------ )   - ------ ( ------ )  ')
    util.disp ('              2 M_AB ( d R_AB )     2 M_BC ( d R_BC )  ')
    util.disp ('                                                       ')
    util.disp ('              cos(theta)    d       d                  ')
    util.disp ('            - ----------  ------  ------               ')
    util.disp ('                 M_B      d R_AB  d R_BC               ')
    util.disp ('                                                       ')
    util.disp ('where M_AB, M_BC are the reduced masses of AB, BC pairs')
    util.disp ('theta is the constant(!) angle between the AB, BC bonds') 
    util.disp ('*******************************************************')
    util.disp ( [ 'Dimensions    : ' num2str(obj.dof)     ] )
    util.disp ( [ 'Mass of atom A: ' num2str(obj.mass(1)) ] )
    util.disp ( [ 'Mass of atom B: ' num2str(obj.mass(2)) ] )
    util.disp ( [ 'Mass of atom C: ' num2str(obj.mass(3)) ] )
    util.disp ( [ 'ABC angle     : ' num2str(obj.theta)   ] )
end

% Check that we have a reasonable grid (i.e. FFT)
if numel(space.dof) ~= 2
	util.error('Triatomic kinetic energy works only for two dimensions.');
end

if ~isa ( space.dof{obj.dof(1)}, 'grid.fft' ) || ...
   ~isa ( space.dof{obj.dof(2)}, 'grid.fft' )
    util.error ('Triatomic kinetic energy only works for fft grids')
end

% Explicitely disable the internal kinetic energy of the grids.
% At least I cannot think of a single case where you would want this
space.dof{obj.dof(1)}.nokin = true;
space.dof{obj.dof(2)}.nokin = true;

% (Reduced) masses
mass_A = obj.mass(1);
mass_B = obj.mass(2);
mass_C = obj.mass(3);

mass_AB = 1 / ( 1/mass_A + 1/mass_B );
mass_BC = 1 / ( 1/mass_B + 1/mass_C );

% Grid representation of kinetic operator. Note d^2/dx^2 => - p_x^2
obj.grid = ...
      space.fbr.grid_ND{obj.dof(1)}.^2 / ( 2 * mass_AB ) ...
    + space.fbr.grid_ND{obj.dof(2)}.^2 / ( 2 * mass_BC ) ...
    + space.fbr.grid_ND{obj.dof(1)}.*space.fbr.grid_ND{obj.dof(2)} ...
    * cos (obj.theta) / mass_B;

obj.grid(obj.grid > hamilt.truncate.delta) = hamilt.truncate.delta;

obj.grid_exp = exp( -1i * time.sub.delta * fraction * obj.grid);
