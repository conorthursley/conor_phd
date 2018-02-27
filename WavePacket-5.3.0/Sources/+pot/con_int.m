% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function con_int
global hamilt space

util.disp (' ')
util.disp ('*************************************************')
util.disp ('Generic E x e conical intersection example :     ')
util.disp ('                                                 ')
util.disp ('Two degenerate electronic states coupled by      ')
util.disp ('two degenerate normal modes of vibration with    ')
util.disp ('linear and quadratic Jahn-Teller coupling        ')
util.disp ('                                                 ')
util.disp ('            omega   ( X^2+Y^2   0  )             ')
util.disp (' V    (R) = ----- * (              )             ')
util.disp ('  dia         2     (  0   X^2+Y^2 )             ')
util.disp ('                                                 ')
util.disp ('                    ( X   Y )                    ')
util.disp ('          + kappa * (       )                    ')
util.disp ('                    ( Y  -X )                    ')
util.disp ('                                                 ')
util.disp ('            gamma   ( +(X^2-Y^2)  -2*X*Y   )     ')
util.disp ('          + ----- * (                      )     ')
util.disp ('              2     (  -2*X*Y   -(X^2-Y^2) )     ')
util.disp ('                                                 ')
util.disp (' The eigenvalues give the 2 adiabatic potentials ')
util.disp ('                                                 ')
util.disp ('                   2       (      2  2           ')
util.disp (' E    (R) = omega R /2 +/- ( kappa  R  + ...     ')
util.disp ('  adi                      (                     ')
util.disp ('                                                 ')
util.disp ('                3                   2  4   ) 1/2 ')
util.disp (' + gamma*kappa*R *cos(3*phi) + gamma  R /4 )     ')
util.disp ('                                           )     ')
util.disp ('                                                 ')
util.disp (' with polar coordinates                          ')
util.disp (' R = ( X^2 + Y^2 )^1/2 and phi = atan(Y/X)       ')
util.disp ('                                                 ')
util.disp (' see e. g. W. Eisfeld and A. Viel                ')
util.disp (' J. Chem. Phys. 122, 204317  (2005) where the    ')
util.disp (' (pseudo) Jahn-Teller is discussed to 6th order !')
util.disp ('                                                 ')
util.disp ('*************************************************')
util.disp ('                                                 ')
util.disp ( [ 'Force constant     omega : ' num2str(hamilt.pot.omega) ] )
util.disp ( [ 'Linear JT-coupling kappa : ' num2str(hamilt.pot.kappa) ] )
util.disp ( [ 'Quadr. JT-coupling gamma : ' num2str(hamilt.pot.gamma) ] )
if space.size.n_dim == 1
    util.disp ( [ 'Energy gap         delta : ' num2str(hamilt.pot.delta) ] )
end
util.disp ( ' ' )

% Check validity
if space.size.n_dim > 2
    util.error ('This potential is only for 1 or 2 dimension')
end

if hamilt.coupling.n_eqs > 2
    util.error ('This potential is only for 1 or 2 equations')
end
  
x = space.dvr.grid_ND{1};
x2 = x.^2;
if space.size.n_dim == 2
    y = space.dvr.grid_ND{2};
    y2 = y.^2;
    r2 = x2 + y2;
else
    r2 = x2;
end

eps1 = hamilt.pot.omega/2 * r2 + hamilt.pot.kappa * x;
eps2 = hamilt.pot.omega/2 * r2 - hamilt.pot.kappa * x;
if space.size.n_dim == 2
    beta  =     + hamilt.pot.kappa   *     y;
    beta = beta - hamilt.pot.gamma   * (x.*y);
    eps1 = eps1 + hamilt.pot.gamma/2 * (x2-y2);
    eps2 = eps2 - hamilt.pot.gamma/2 * (x2-y2);
else
    beta  =     + hamilt.pot.gamma;
    eps1 = eps1 + hamilt.pot.delta/2;
    eps2 = eps2 - hamilt.pot.delta/2;
end

% Coupled two state problem
if hamilt.coupling.n_eqs ==2
    hamilt.pot.grid_ND{1,1} = eps1;
    hamilt.pot.grid_ND{2,2} = eps2;
    hamilt.pot.grid_ND{1,2} = beta;
    
% Only lower adiabatic state: Mexican hat
elseif hamilt.coupling.n_eqs ==1
    eta   = (eps1+eps2)/2;
    delta = (eps1-eps2)/2;
    rho   = sqrt(delta.^2+beta.^2);
    hamilt.pot.grid_ND{1,1} = eta - rho;
end
