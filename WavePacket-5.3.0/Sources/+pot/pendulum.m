% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015 Burkhard Schmidt
%
% see the README file for license details.

function pendulum
global hamilt space

util.disp (' ')
util.disp ('************************************************************')
util.disp ('Potential for generalized quantum pendula                   ')
util.disp ('                                                            ')
util.disp (' V = xi*csc^2 + chi*cot*csc - eta*cos - zeta*cos^2 + v_0    ')
util.disp (' see work by B. Schmidt and B. Friedrich (2014/2015/2017)   ')
util.disp ('                                                            ')
util.disp ('Important note on usage                                     ')
util.disp ('Planar pendulum:    0 <= theta <= 2*pi (or also 4*pi !)     ')
util.disp ('Spherical pendulum: 0 <= theta <= pi                        ')
util.disp ('                    and xi = m^2 - 1/4  where m is the      ')
util.disp ('                    azimuthal quantum number                ')
util.disp ('************************************************************')

if ~isfield (hamilt.pot,'xi')   hamilt.pot.xi=0;   end
if ~isfield (hamilt.pot,'chi')  hamilt.pot.chi=0;  end
if ~isfield (hamilt.pot,'eta')  hamilt.pot.eta=0;  end
if ~isfield (hamilt.pot,'zeta') hamilt.pot.zeta=0; end
if ~isfield (hamilt.pot,'v_0')  hamilt.pot.v_0=0;  end
util.disp ( [ 'csc^2 term  (xi)   : ' num2str(hamilt.pot.xi) ] )
util.disp ( [ 'cot*csc term (chi) : ' num2str(hamilt.pot.chi) ] )
util.disp ( [ 'Orientation (eta)  : ' num2str(hamilt.pot.eta) ] )
util.disp ( [ 'Alignment (zeta)   : ' num2str(hamilt.pot.zeta) ] )
util.disp ( [ 'Energy shift (v_0) : ' num2str(hamilt.pot.v_0) ] )

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential only for one dimension')
end

if hamilt.coupling.n_eqs ~= 1
    util.error ('This potential only for single Schrödinger equation')
end

% Truncate csc^2 = 1/sin^2 to prevent singularity
my_eps = 1e-5;
x = space.dvr.grid_ND{1};
csc2 = zeros(size(space.dvr.grid_ND{1}));
for ii=1:length(x)
    if abs(x(ii)) < my_eps || abs(x(ii)-pi) < my_eps
        csc2(ii) = 1/sin(my_eps)^2;
    else
        csc2(ii) = 1/sin(x(ii))^2;
    end
end
        
% Evaluate potential energy function
hamilt.pot.grid_ND{1,1} = ...
    + hamilt.pot.xi*csc2 ...
    + hamilt.pot.chi*csc2.*cos(x) ...
    - hamilt.pot.eta*cos(x) ...
    - hamilt.pot.zeta*cos(x).^2 ...
    + hamilt.pot.v_0;









