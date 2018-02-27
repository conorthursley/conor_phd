% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function tully3 
global hamilt space

util.disp (' ')
util.disp ('******************************************************************')
util.disp ('Extended coupling example: J.C.Tully, J.Chem.Phys. 93, 1061 (1990)')
util.disp ('******************************************************************')
util.disp ( [ 'Tully parameter A : ' num2str(hamilt.pot.A) ] )
util.disp ( [ 'Tully parameter B : ' num2str(hamilt.pot.B) ] )
util.disp ( [ 'Tully parameter C : ' num2str(hamilt.pot.C) ] )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential is only for 1 dimension')
end

if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential is only for 2 states')
end
    
hamilt.pot.grid_ND{1,1} = + hamilt.pot.A * ones(size(space.dvr.grid_ND{1}));
hamilt.pot.grid_ND{2,2} = - hamilt.pot.A * ones(size(space.dvr.grid_ND{1}));

left  = find ( space.dvr.grid_ND{1}< 0 );
right = find ( space.dvr.grid_ND{1}>=0 );
hamilt.pot.grid_ND{1,2}(left)  = + hamilt.pot.B *       exp (   hamilt.pot.C * space.dvr.grid_ND{1}(left)  );
hamilt.pot.grid_ND{1,2}(right) = + hamilt.pot.B * ( 2 - exp ( - hamilt.pot.C * space.dvr.grid_ND{1}(right) ) );
hamilt.pot.grid_ND{1,2} = hamilt.pot.grid_ND{1,2}';


