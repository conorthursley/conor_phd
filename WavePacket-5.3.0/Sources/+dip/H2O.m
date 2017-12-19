% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Ulf Lorenz
%
% see the README file for license details.

function H2O
global hamilt space

util.disp (' ')
util.disp ('*****************************************************')
util.disp ('Model function for transition dipole moment of H2O.  ')
util.disp ('See J. Chem. Phys. 97:8285                           ')
util.disp ('*****************************************************')

%% Check validity
if space.size.n_dim ~= 2
    util.error ('This dipole function is only for 2 dimensions')
end

if hamilt.coupling.n_eqs ~= 2
    util.error ('This dipole function is only for 2 states')
end
 
%% Set parameters
b  = 1;     % range parameter
r0 = 1.81;  % equilibrium distance H-O

% dipole moment along x, so use a laser with polar = 0.
hamilt.d_x.grid_ND{1,2} = 2.225^2 ./ ( ...
        (1 + exp(b * (space.dvr.grid_ND{1} - r0))) .* ...
        (1 + exp(b * (space.dvr.grid_ND{2} - r0))) );
    
hamilt.d_x.grid_ND{1,1} = [];
hamilt.d_x.grid_ND{2,2} = [];


% No dipole moments along y
hamilt.d_y.grid_ND{1,1} = [];
hamilt.d_y.grid_ND{2,2} = [];
hamilt.d_y.grid_ND{1,2} = [];


