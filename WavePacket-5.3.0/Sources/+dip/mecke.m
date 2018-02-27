% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function mecke
global hamilt space

util.disp (' ')
util.disp ('*****************************************************')
util.disp ('Mecke model function for permanent dipole moment     ')
util.disp ('                                                     ')
util.disp ('  mu(r) = q_0 * r * exp ( - r / r_0 )                ')
util.disp ('                                                     ')
util.disp ('   see, e.g., Jakubetz, Manz, Mohan                  ')
util.disp ('   J. Chem. Phys. 90, 3686 (1989)                    ')
util.disp ('*****************************************************')
util.disp ( [ 'Charge parameter q_0: ' num2str(hamilt.dip.q_0) ] )
util.disp ( [ 'Length parameter r_0: ' num2str(hamilt.dip.r_0) ] )
util.disp ( ' ' )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This dipole function is only for 1 dimension')
end

if hamilt.coupling.n_eqs ~= 1
    util.error ('This dipole function is only for 1 state')
end
    
% Mecke function (dipole moment along x)
hamilt.d_x.grid_ND{1,1} = hamilt.dip.q_0 * space.dvr.grid_ND{1}...
    .* exp ( - space.dvr.grid_ND{1}/hamilt.dip.r_0 );

% No dipole moment along y
hamilt.d_y.grid_ND{1,1} = [];
