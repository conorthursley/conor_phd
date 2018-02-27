% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function metiu1

global hamilt space

util.disp(' ')
util.disp('**********************************************************')
util.disp('Model potential: Dateo&Metiu, J.Chem.Phys. 95, 7392 (1991)')
util.disp(' ')
util.disp(' ')
util.disp(' V(\Theta) = A  + A P (cos \Theta)  + A P (cos \Theta)    ')
util.disp('              0    1 1                 2 2                ')
util.disp('**********************************************************')
util.disp([ 'Parameter A0 : ' num2str(hamilt.pot.a0) ])
util.disp([ 'Parameter A1 : ' num2str(hamilt.pot.a1) ])
util.disp([ 'Parameter A2 : ' num2str(hamilt.pot.a2) ])

if space.size.n_dim > 1 
    util.err('This potential works only for one DOF')
end

if ~isa(space.dof{1}, 'grid.legendre')
    util.err('This potential works only for Legendre DVR')
end

if hamilt.coupling.n_eqs > 1
    util.err('This potential works only for one Schroedinger equation')
end

% Note that a call to legendre(N,x) returns all associated Legendre polynomials
% of degree N. So we have to mask out the values for m=0 be hand and recast them
% to the original form, which is a bit of annoying.
leg1 = legendre(1, space.dvr.grid_ND{1}(:));
leg1 = leg1(1, :)';
leg2 = legendre(2, space.dvr.grid_ND{1}(:));
leg2 = leg2(1, :)';

hamilt.pot.grid_ND{1,1} = hamilt.pot.a0 + hamilt.pot.a1 * leg1 ...
                          + hamilt.pot.a2 * leg2;
