%--------------------------------------------------------------------------
%
% Set up grid representation of p-th additional multiplicative operators
% Rotational characteristics using cos^n as projector
% 
% n=1: Orientation
% n=2: Alignment
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2008 Ulf Lorenz
%               2008-2016 Burkhard Schmidt
%
% see the README file for license details.

function cosine (p)

global hamilt space

util.disp ( '****************************************************' )
util.disp ([ int2str(p) '-th additional multiplicative operators:' ] )
util.disp ( 'Directional property (cos^n) in one dimension' )
util.disp ( '****************************************************' )

% Set/check/print input parameters
if ~isfield(space.amo{p},'exp')
    space.amo{p}.exp = 1;
end

if ~isfield(space.amo{p},'dof')
    space.amo{p}.dof = 1;
end

if ~isfield(space.amo{p},'label')
    switch space.amo{p}.exp
        case 1
            space.amo{p}.label = 'Orientation';
        case 2
            space.amo{p}.label = 'Alignment';
        otherwise
            space.amo{p}.label = ['<cos^' int2str(space.amo{p}.exp) '\theta>'];
    end
end

% Output
util.disp ( [ 'Exponent n       : ' num2str(space.amo{p}.exp)])
if space.size.n_dim > 1
    util.disp ( [ 'Degree of freedom: ' space.amo{p}.dof] )
end
util.disp ( [ 'Label            : ' space.amo{p}.label])

% Use same projector for all (coupled) states
for m=1:hamilt.coupling.n_eqs
    
    % Create the projection grid
    if isa(space.dof{space.amo{p}.dof}, 'grid.legendre')
        util.disp( 'Treating position variable as cos \theta' )
        space.amo{p}.grid_ND{m,m} =     space.dvr.grid_ND{space.amo{p}.dof} .^space.amo{p}.exp;
    else
        space.amo{p}.grid_ND{m,m} = cos(space.dvr.grid_ND{space.amo{p}.dof}).^space.amo{p}.exp;
    end
    
end
util.disp (' ')

