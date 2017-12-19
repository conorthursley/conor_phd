%--------------------------------------------------------------------------
%
% Set up grid representation of p-th additional multiplicative operators
% Gaussian bell-shaped function as projector 
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function gauss (p)

global hamilt space

util.disp ( '****************************************************' )
util.disp ([ int2str(p) '-th additional multiplicative operators:' ] )
util.disp ( 'Gaussian bell-shaped projection function in N-dim' )
util.disp ( '                                                    ' )
util.disp ( '           N      [   ( Ri-R0i )^2 ]                ' )
util.disp ( ' G (R) = Prod exp [ - (--------)   ]                ' )
util.disp ( '          i=1     [   (  2*Wi  )   ]                ' )
util.disp ( '                                                    ' )
util.disp ( 'where the product extends over all dimensions       ' )
util.disp ( '                                                    ' )
util.disp ( 'For an example, see Eq. (49) in JCP 109, 385 (1998) ' )
util.disp ( 'DOI:10.1063/1.476575 by W. Zhu and H. Rabitz        ' )
util.disp ( '****************************************************' )

% Check input parameters
if length(space.amo{p}.pos_0)~=space.size.n_dim
    util.error ('Incompatible dimensionality for pos_0 of projection')
end

if length(space.amo{p}.width)~=space.size.n_dim
    util.error ('Incompatible dimensionality for width of projection')
end

if ~isfield(space.amo{p},'label')
    space.amo{p}.label = 'Gaussian';
end

% Output
util.disp ( [ 'Mean value position       R0 : ' num2str(space.amo{p}.pos_0) ] )
util.disp ( [ 'Position uncertainty      W  : ' num2str(space.amo{p}.width) ] )
util.disp ( [ 'Label                        : '         space.amo{p}.label  ] )
util.disp (   ' ')

% Use same projector for all (coupled) states
for m=1:hamilt.coupling.n_eqs
    space.amo{p}.grid_ND{m,m} = ones ( size(space.dvr.grid_ND{1}) );
    
    % Tensor product of one-dimensional Gaussians
    for k = 1:space.size.n_dim
        space.amo{p}.grid_ND{m,m} = space.amo{p}.grid_ND{m,m} / ...
		(2*space.amo{p}.width(k)*sqrt(pi)) .* ...
		exp (  -((space.dvr.grid_ND{k}-space.amo{p}.pos_0(k)) / (space.amo{p}.width(k)*2)).^2  );
    end
    
end


