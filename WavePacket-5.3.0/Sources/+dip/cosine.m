% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.


function cosine
global hamilt space

util.disp (' ')
util.disp ('*****************************************************')
util.disp ('Permanent dipole moment as a power of cos Theta      ')
util.disp (' mu(Theta) = f * cos^n  Theta                        ')
util.disp ('*****************************************************')

% Set default values
if ~isfield(hamilt.dip, 'pre')
        hamilt.dip.pre = 1;
end
util.disp (['Prefactor f : ' num2str(hamilt.dip.pre)])

if ~isfield(hamilt.dip, 'exp')
        hamilt.dip.exp = 1;
end
util.disp (['Exponent  n : ' num2str(hamilt.dip.exp)])

% Check/print input parameters
if hamilt.coupling.n_eqs ~= 1
    util.error ('This dipole function is only for 1 state')
end
    
if space.size.n_dim > 1
    util.error ('This dipole function is only for 1 dof')
end

% Set up the grid representation of dipole moment along x
if isa (space.dof{space.amo{1}.dof}, 'grid.legendre')
    util.disp( 'Treating position variable as cos \theta' )
    hamilt.d_x.grid_ND{1,1} = hamilt.dip.pre *     space.dvr.grid_ND{1} .^hamilt.dip.exp;
else
    hamilt.d_x.grid_ND{1,1} = hamilt.dip.pre * cos(space.dvr.grid_ND{1}).^hamilt.dip.exp;
end

% No dipole moment along y
hamilt.d_y.grid_ND{1,1} = [];
