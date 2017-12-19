% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015 Burkhard Schmidt
%
% see the README file for license details.

function razavy

global hamilt space

util.disp (' ')
util.disp ('*******************************************************')
util.disp ('Potential energy for a Razayv single/double well       ')
util.disp ('American Journal of Physics, 48(4), 285-288 (1980)     ')
util.disp ('DOI:10.1119/1.12141, Eq. (2.2) with beta=1             ')
util.disp ('                                                       ')
util.disp (' V (R) = - s*(n+1)*cosh(2*x) + s^2/8*cosh(4*x) - s^2/8 ')
util.disp ('                                                       ')
util.disp (' or our modified version (Mirahmadi/Schmidt 2016/17)   ')
util.disp ('                                                       ')
util.disp (' V (R) = eta*cosh(x) + zeta*cosh^2(x)                  ')
util.disp ('                                                       ')
util.disp ('*******************************************************')

% Check validity
if space.size.n_dim ~= 1
    util.error ('This potential only for one dimension')
end
x = space.dvr.grid_ND{1};

if hamilt.coupling.n_eqs ~= 1
    util.error ('This potential only for single Schrödinger equation')
end

if ~isfield (hamilt.pot,'modified')
    hamilt.pot.modified = false;
end

if ~hamilt.pot.modified
    
    s = hamilt.pot.s;
    n = hamilt.pot.n;
    
    util.disp ('Original version of Razayv single/double well potential')
    util.disp ( [ 'Strength parameter  s : ' num2str(s) ] )
    util.disp ( [ 'Magic number n        : ' num2str(n) ] )
    
    hamilt.pot.grid_ND{1,1} = ...
        - s*(n+1)*cosh(2*x) ...
        + s^2/8*cosh(4*x) ...
        - s^2/8;
    
else
    
    eta   = hamilt.pot.eta;
    zeta  = hamilt.pot.zeta;
    beta  = sqrt(zeta)*sign(eta);
    kappa = abs(eta/beta);
    
    util.disp ('Modified version of Razayv single/double well potential')
    util.disp ( [ 'Prefactor of cosh    (eta) : ' num2str(eta) ] )
    util.disp ( [ 'Prefactor of cosh^2 (zeta) : ' num2str(zeta) ] )
    util.disp ( [ 'SUSY parameter      (beta) : ' num2str(beta) ] )
    util.disp ( [ 'Topological index  (kappa) : ' num2str(kappa) ] )
    util.disp ( [ 'Energy at minima           : ' num2str(-eta^2/(4*zeta)) ] )
    util.disp ( [ 'Position of minima         : ' num2str(acosh(-eta/(2*zeta))) ] )
    util.disp ( [ 'Energy at maximum          : ' num2str(eta+zeta) ] )
    
    hamilt.pot.grid_ND{1,1} = eta * cosh(x) + zeta * cosh(x).^2;
    
end






