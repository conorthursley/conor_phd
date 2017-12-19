%--------------------------------------------------------------------------
%
% Set up grid representation of p-th additional multiplicative operators
% Educts=Reactant versus Products in a chemical exchange reaction
%  A + BC -> AB + C       
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%               2008-2016 Burkhard Schmidt
%
% see the README file for license details.


function reaction (p)

global hamilt space

util.disp ( '****************************************************')
util.disp ([ int2str(p) '-th additional multiplicative operators:' ] )
util.disp ('Projection on reactant/product region   ')
util.disp ('for a chemical exchange reaction     ')
util.disp ('          A + BC -> AB + C           ')
util.disp ('****************************************************')

% Set/check/print input parameters
if space.size.n_dim < 2
    util.error('At least two degrees of freedom required!');
end

if ~isfield(space.amo{p},'reac')
    space.amo{p}.reac = 1;
end
util.disp(['Index of reactant distance coordinate AB :' num2str(space.amo{p}.reac)])

if ~isfield(space.amo{p},'prod')
    space.amo{p}.prod = 2;
end
util.disp(['Index of product distance coordinate BC :' num2str(space.amo{p}.prod)])

if ~isfield(space.amo{p},'side')
    space.amo{p}.side = 'p';
end
switch lower(space.amo{p}.side)
    case 'r'
        util.disp('Projecting on reactant part')
    case 'p'
        util.disp('Projecting on product part')
    otherwise
        util.error('space.amo{p}.side has to be "r" or "p"')
end
util.disp( ' ')

if ~isfield(space.amo{p},'label')
    switch lower(space.amo{p}.side)
    case 'r'
        space.amo{p}.label = 'Reactant';
    case 'p'
        space.amo{p}.label = 'Product';
    end
end
util.disp ( [ 'Label                            : '         space.amo{p}.label] )

% Use same projector for all (coupled) states
for m=1:hamilt.coupling.n_eqs
    
    switch lower(space.amo{p}.side)
        case 'r'
            space.amo{p}.grid_ND{m,m} = ...
                (space.dvr.grid_ND{space.amo{p}.reac} > space.dvr.grid_ND{space.amo{p}.prod});
        case 'p'
            space.amo{p}.grid_ND{m,m} = ...
                (space.dvr.grid_ND{space.amo{p}.reac} < space.dvr.grid_ND{space.amo{p}.prod});
    end
    
end
