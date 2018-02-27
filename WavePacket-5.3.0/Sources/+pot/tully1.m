% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function tully1
global hamilt space

util.disp (' ')
util.disp ('****************************************************************')
util.disp ('Single crossing example: J.C.Tully, J.Chem.Phys. 93, 1061 (1990)')
util.disp ('****************************************************************')
util.disp ( [ 'Tully parameter A : ' num2str(hamilt.pot.A) ] )
util.disp ( [ 'Tully parameter B : ' num2str(hamilt.pot.B) ] )
util.disp ( [ 'Tully parameter C : ' num2str(hamilt.pot.C) ] )
util.disp ( [ 'Tully parameter D : ' num2str(hamilt.pot.D) ] )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential is only for 1 dimension')
end

if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential is only for 2 states')
end
    
hamilt.pot.grid_ND{1,1} = hamilt.pot.A * sign (space.dvr.grid_ND{1}) .* ( 1 - exp(-hamilt.pot.B*sign (space.dvr.grid_ND{1}).*space.dvr.grid_ND{1}) );
hamilt.pot.grid_ND{2,2} = - hamilt.pot.grid_ND{1,1};
hamilt.pot.grid_ND{1,2} = hamilt.pot.C * exp ( - hamilt.pot.D * space.dvr.grid_ND{1} .^ 2);



