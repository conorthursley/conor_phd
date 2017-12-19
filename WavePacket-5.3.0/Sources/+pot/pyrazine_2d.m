% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function pyrazine_2d
global hamilt space

util.disp (' ')
util.disp ('*************************************************')
util.disp ('Ultrafast S2-S1 internal conversion of pyrazine  ')
util.disp ('Two state / two mode model exhibiting conical    ')
util.disp ('intersections along a straight line in the       ')
util.disp ('(R1,R2) plane. Note that these intersections     ')
util.disp ('are not(!) due to a Jahn-Teller effect but are   ')
util.disp ('only partially symmetry induced!                 ')
util.disp ('                                                 ')
util.disp (' Reduced to two dimensions setting R3=const      ')
util.disp ('                                                 ')
util.disp ('Two non-degenerate electronic states (D_2h group)')
util.disp ('    S1 (B3u)                                     ')
util.disp ('    S2 (B2u)                                     ')
util.disp ('Two tuning modes (R1,R2), one coupling mode (R3) ')
util.disp ('    R1 = mode 01  (totally symmetric)            ') 
util.disp ('    R2 = mode 06a (totally symmetric)            ')
util.disp ('    R3 = mode 10a = const (B1g symmetry)         ') 
util.disp ('                                                 ')
util.disp ('             3  omega_k   ( R_k^2   0  )         ')
util.disp (' V    (R) = Sum ------- * (            )         ')
util.disp ('  dia       k=1   2       (  0   R_k^2 )         ')
util.disp ('                                                 ')
util.disp ('             2  ( kappa1_i R_i      0        )   ')
util.disp ('          + Sum (                            )   ')
util.disp ('            i=1 (     0         kappa2_i R_i )   ')
util.disp ('                                                 ')
util.disp ('             3             (  0   R_j )          ')
util.disp ('          + Sum lambda12 * (          )          ')
util.disp ('            j=3            ( R_j   0  )          ')
util.disp ('                                                 ')
util.disp ('            ( E1   0  )                          ')
util.disp ('          + (         )                          ')
util.disp ('            ( 0    E2 )                          ')
util.disp ('                                                 ')
util.disp (' Note that the summations over kappa and lambda  ')
util.disp (' extend over tuning and coupling modes, respect.,')
util.disp (' while the omega summation extends over all modes')
util.disp ('                                                 ')
util.disp (' Schneider, Domcke,         CPL 150, 235 (1988)  ')
util.disp (' Schneider, Domcke, Köppel, JCP 92, 1045 (1990)  ')
util.disp ('                                                 ')
util.disp ('*************************************************')
util.disp ('                                                 ')
util.disp ( [ 'Force constant omega (S0=S1=S2): ' num2str(hamilt.pot.omega ) ] )
util.disp ( [ 'Linear parameters kappa (S1)   : ' num2str(hamilt.pot.kappa1) ] )
util.disp ( [ 'Linear parameters kappa (S2)   : ' num2str(hamilt.pot.kappa2) ] )
util.disp ( [ 'Linear coupling lambda (S1-S2) : ' num2str(hamilt.pot.lambda) ] )
util.disp ( [ 'Vertical excitation E (S1)     : ' num2str(hamilt.pot.energ1) ] )
util.disp ( [ 'Vertical excitation E (S2)     : ' num2str(hamilt.pot.energ2) ] )
util.disp ( [ 'Fixed value of coupling mode R3: ' num2str(hamilt.pot.R3   ) ] )
util.disp ( ' ' )

% Check validity
if space.size.n_dim ~= 2
    util.error ('This potential is only for 2 dimension')
end

if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential is only for 2 equations')
end
  

hamilt.pot.grid_ND{1,1} = hamilt.pot.omega(1)/2 * space.dvr.grid_ND{1}.^2 ...
                     + hamilt.pot.omega(2)/2 * space.dvr.grid_ND{2}.^2 ...
                     + hamilt.pot.kappa1(1)  * space.dvr.grid_ND{1} ...
                     + hamilt.pot.kappa1(2)  * space.dvr.grid_ND{2} ...
                     + hamilt.pot.energ1;

hamilt.pot.grid_ND{2,2} = hamilt.pot.omega(1)/2 * space.dvr.grid_ND{1}.^2 ...
                     + hamilt.pot.omega(2)/2 * space.dvr.grid_ND{2}.^2 ...
                     + hamilt.pot.kappa2(1)  * space.dvr.grid_ND{1} ...
                     + hamilt.pot.kappa2(2)  * space.dvr.grid_ND{2} ...
                     + hamilt.pot.energ2;
                 
hamilt.pot.grid_ND{1,2} = hamilt.pot.lambda     * hamilt.pot.R3;
