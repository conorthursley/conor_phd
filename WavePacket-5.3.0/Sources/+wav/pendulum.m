%------------------------------------------------------------------------------
%
% Ordinary (angular) Mathieu function to be used as initial wavefunction 
% for quantum problems defined on a ring
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008 Ulf Lorenz
%               2008-2009 Burkhard Schmidt
%
% see the README file for license details.

function init_grid = pendulum(dir)

global psi space

util.disp (' ')
util.disp ('******************************************************' )
util.disp ('Initial wavefunction of a plane pendulum:             ' )
util.disp ('                                                      ' )
util.disp ('Ordinary (angular) Mathieu functions                  ' )
util.disp ('                                                      ' )
util.disp (' ce_r (eta,q) with r>=0 : cosine elliptic functions   ' )
util.disp (' se_r (eta,q) with r>=1 :   sine elliptic functions   ' )
util.disp ('                                                      ' )
util.disp ('solving a Schrödinger equation for potential          ' )
util.disp ('                                                      ' )
util.disp ('     V(phi) = V0/2 [ 1 + cos (m*(R-R0) ]              ' )
util.disp ('                                                      ' )
util.disp ('for periodic coordinate 0 < R <= 2*pi where the       ' )
util.disp ('multiplicity m indicates the number of maxima/minima  ' )
util.disp ('and where eta=m*(R-R0)/2 and q=2V0/m^2                ' )
util.disp ('                                                      ' )
util.disp ('m=1: Only pi-periodic solutions: Only even orders r   ' )
util.disp ('m=2: All 2pi-periodic solutions: Only int. orders r   ' )
util.disp ('m>2: Requires also fractional orders of r             ' )
util.disp ('     Yields complex solutions                         ' )
util.disp ('******************************************************' )

init_grid = ones (size(space.dvr.grid_ND{1}));

% Check for correct periodicity
switch psi.dof{dir}.multiple
    case 1 % Only even integer orders r
        if psi.dof{dir}.order~=2*round(psi.dof{dir}.order/2)
            util.error ('Order of Mathieu functions not compatible with m!')
        end
    case 2 % Only integer orders r 
        if psi.dof{dir}.order~=round(psi.dof{dir}.order)
            util.error ('Order of Mathieu functions not compatible with m!')
        end            
    otherwise
        util.error ('Fractional order Mathieu functions not yet implemented!')
end
     
% Arguments of Mathieu equation
r = psi.dof{dir}.order;
q = psi.dof{dir}.barrier*2/psi.dof{dir}.multiple^2;
e = psi.dof{dir}.multiple*(space.dvr.grid_ND{dir}-psi.dof{dir}.shift)/2;

switch lower(psi.dof{dir}.parity)
    case ('c')
        util.disp ( 'Cosine elliptic' )
        init_grid = init_grid .* math.cee(r,q,e);
    case {'s'}
        util.disp ( 'Sine elliptic' )
        init_grid = init_grid .* math.see(r,q,e);
    case {'l'}
        util.disp ( 'Localized combination of cosine and sine elliptic' )
        init_grid = init_grid .* (math.cee(r,q,e)+math.see(r+1,q,e))/sqrt(2);
    otherwise
        util.error ('Wrong choice for parity of Mathieu function')
end
util.disp ( [ 'Mathieu order r      : ' int2str( psi.dof{dir}.order) ] )
util.disp ( [ 'Multiplicity m       : ' int2str( psi.dof{dir}.multiple) ] )
util.disp ( [ 'Potential barrier V0 : ' num2str( psi.dof{dir}.barrier) ] )
util.disp ( [ 'Potential shift R0   : ' num2str( psi.dof{dir}.shift) ] )




