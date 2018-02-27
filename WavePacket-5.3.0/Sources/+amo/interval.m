%--------------------------------------------------------------------------
%
% Set up grid representation of p-th additional multiplicative operators
% Characteristic function on an interval (in 1 dimension)
% Characteristic function on a rectangle (in 2 dimensions)
% etc ...
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function interval (p)

global hamilt space

util.disp ( '****************************************************' )
util.disp ([ int2str(p) '-th additional multiplicative operators:' ] )
util.disp ( 'Interval / rectangle / cuboid / ...' )
util.disp ( '****************************************************' )

% Check input parameters
if length(space.amo{p}.min)~=space.size.n_dim
    util.error ('Incompatible dimensionality for beginning of projection')
end

if length(space.amo{p}.max)~=space.size.n_dim
    util.error ('Incompatible dimensionality for end of projection')
end

if any (space.amo{p}.max <= space.amo{p}.min)
    util.error ( 'Wrong ordering of projection interval min/max parameters' )
end

if ~isfield(space.amo{p},'label')
    space.amo{p}.label = 'Interval';
end

% Output
util.disp ( [ 'Beginning of projection interval : ' num2str(space.amo{p}.min) ] )
util.disp ( [ 'End of projection interval       : ' num2str(space.amo{p}.max) ] )
util.disp ( [ 'Label                            : '         space.amo{p}.label] )
util.disp (   ' ')

% Use same projector for all (coupled) states
for m=1:hamilt.coupling.n_eqs
    space.amo{p}.grid_ND{m,m} = ones(size(space.dvr.grid_ND{1}));
    
    % Tensor product of one-dimensional Gaussians
    for k = 1:space.size.n_dim
        
        % Find grid points outside interval and set AMO function to zero
        space.amo{p}.grid_ND{m,m}(space.dvr.grid_ND{k} < space.amo{p}.min(k) ...
            | space.dvr.grid_ND{k} > space.amo{p}.max(k)) = 0;
        
    end
    
end
