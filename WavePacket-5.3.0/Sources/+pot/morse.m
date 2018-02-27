% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2009 Burkhard Schmidt
%
% see the README file for license details.

function morse
global hamilt space 

util.disp (' ')
util.disp ('*******************************************************')
util.disp ('Potential energy: Morse oscillator')
util.disp ('*******************************************************')
for k=1:space.size.n_dim
    util.disp ( [ 'Dimension : ' num2str(k) ])
    util.disp ( [ 'Dissociation energy    : ' num2str(hamilt.pot.d_e(k,:)) ] )
    util.disp ( [ 'Equilibrium position   : ' num2str(hamilt.pot.r_e(k,:)) ] )
    util.disp ( [ 'Range parameter (alfa) : ' num2str(hamilt.pot.alf(k,:)) ] )
    util.disp ('')
end
if ~isfield (hamilt.pot,'t_e')
    hamilt.pot.t_e = 0;
end
util.disp ( [ 'Energetic shifts    : ' num2str(hamilt.pot.t_e) ] )
util.disp (' ')

% Loop over coupled equations
for m=1:hamilt.coupling.n_eqs

    % Summing up contributions from each component of position vector
    hamilt.pot.grid_ND{m,m} = zeros(size(space.dvr.grid_ND{1}));

    for k = 1:space.size.n_dim
        hamilt.pot.grid_ND{m,m} = hamilt.pot.grid_ND{m,m} + ...
            hamilt.pot.d_e(k,m) * ( 1 - exp ( -hamilt.pot.alf(k,m) * ( space.dvr.grid_ND{k}-hamilt.pot.r_e(k,m) ) ) ).^2 + hamilt.pot.t_e(m);
    end


end

