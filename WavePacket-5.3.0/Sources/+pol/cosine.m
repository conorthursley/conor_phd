% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.


function cosine
global hamilt space

util.disp (' ')
util.disp ('**************************************')
util.disp ('Polarizability as a power of cos Theta')
util.disp (' alpha(Theta) = f * cos^n Theta       ')
util.disp ('**************************************')

% Set default values
if ~isfield(hamilt.pol, 'pre')
        hamilt.pol.pre = 1;
end
util.disp (['Prefactor f : ' num2str(hamilt.pol.pre)])

if ~isfield(hamilt.pol, 'exp')
        hamilt.pol.exp = 2;
end
util.disp (['Exponent  n : ' num2str(hamilt.pol.exp)])

% Check/print input parameters
if hamilt.coupling.n_eqs ~= 1
    util.error ('This polarizability is only for 1 state')
end

if space.size.n_dim > 1
    util.error ('This polarizability is only for 1 dof')
end

% Set up the grid representation along x
if isa (space.dof{space.amo{1}.dof}, 'grid.legendre')
    util.disp( 'Treating position variable as cos \theta' )
    hamilt.p_x.grid_ND{1,1} = hamilt.pol.pre *     space.dvr.grid_ND{1} .^hamilt.pol.exp;
else
    hamilt.p_x.grid_ND{1,1} = hamilt.pol.pre * cos(space.dvr.grid_ND{1}).^hamilt.pol.exp;
end

% No polarizability along y
 hamilt.p_y.grid_ND{1,1} = [];