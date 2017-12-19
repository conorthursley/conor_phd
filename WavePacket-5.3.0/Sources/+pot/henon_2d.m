% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function henon_2d
global hamilt space

util.disp (' ')
util.disp ('*******************************************************')
util.disp ('Potential energy: Henon-Heiles system in 2 dimensions  ')
util.disp ('                                                       ')
util.disp ('         A    2    2         2       3                 ')
util.disp (' V (R) = - ( R  + R ) + L ( R  R  - R /3 )             ')
util.disp ('         2    1    2         1  2    2                 ')
util.disp ('                                                       ')
util.disp (' or in polar coordinates                               ')
util.disp ('                                                       ')
util.disp ('        A    2   1          3                          ')
util.disp (' V(R) = - |R|  + - L |R| sin (3 theta)                 ')
util.disp ('        2        3                                     ')
util.disp ('                                                       ')
util.disp (' which makes the C_3v symmetry obvious. In addition to ')
util.disp (' the minimum found at |R|=0, there are three saddles   ')
util.disp (' at |R| = A/L and theta = pi/3, pi, 5*pi/3 with        ')
util.disp (' a corresponding energy of A^3/(6*L^2).                ')
util.disp (' Note that a classical particle can escape to infinity ')
util.disp (' if its energy exceeds that of the three saddle points.')
util.disp (' Quantum mechanical bound states can only be found     ')
util.disp (' below the energy of the saddles. Strictly speaking,   ')
util.disp (' these states are resonances because they can decay    ')
util.disp (' by tunneling through the barriers.                    ')
util.disp ('                                                       ')
util.disp (' M.J.Davis, E.J.Heller, J. Chem. Phys. 71, 3383 (1979) ')
util.disp (' J.-P. Kuska / C. Herrmann, Physik Journal, Sept. 2005 ')
util.disp ('                                                       ')
util.disp ('*******************************************************')
util.disp ( [ 'Force constant A       : ' num2str(hamilt.pot.A) ] )
util.disp ( [ 'Cubic parameter lambda : ' num2str(hamilt.pot.L) ] )
util.disp ( [ 'Saddle height          : ' num2str(hamilt.pot.A^3/(6*hamilt.pot.L^2)) ] )
util.disp ( [ 'Radius of triangle     : ' num2str(hamilt.pot.A  /   hamilt.pot.L   ) ] )

% Check validity
if hamilt.coupling.n_eqs ~= 1
    util.error ('This potential only for single Schrödinger equation')
end

if space.size.n_dim ~= 2
    util.error ('This potential only for two dimensions')
end

hamilt.pot.grid_ND{1,1} = ...
    + hamilt.pot.A/2 ...
    * ( space.dvr.grid_ND{1}.^2 ...
    + space.dvr.grid_ND{2}.^2 ) ...
    + hamilt.pot.L ...
    * ( space.dvr.grid_ND{1}.^2 .* space.dvr.grid_ND{2}  ...
    - space.dvr.grid_ND{2}.^3  / 3);








