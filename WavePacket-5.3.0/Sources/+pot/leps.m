% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

function leps

global space hamilt

util.disp( '')
util.disp( '********************************************************' )
util.disp( 'LEPS potential for atom/diatom A + BC reactive scattering ' )
util.disp( 'LEPS = London, Eyring, Polanyi, Sato                    ' )
util.disp( 'First two bond coordinates are assumed to be coordinates' )
util.disp( 'corresponding to the two degrees of freedom             ' )
util.disp( '********************************************************' )

util.disp( ['Bending angle       : ' num2str(hamilt.pot.theta)] )
for ii = 1:3
    util.disp( '' )
    util.disp( ['Bond length #       : ' num2str(ii)                       ] )
    util.disp( ['dissociation energy : ' num2str(hamilt.pot.d_e(ii))] )
    util.disp( ['equilibrium distance: ' num2str(hamilt.pot.r_e(ii))] )
    util.disp( ['range parameter     : ' num2str(hamilt.pot.a(ii))  ] )
    util.disp( ['Sato parameter      : ' num2str(hamilt.pot.s(ii))  ] )
end

%% First some consistency checks

if hamilt.coupling.n_eqs > 1
    util.error( 'This surface is only for one Schroedinger equation.' )
end

if length(space.dof) ~= 2
    util.error( 'This surface is for 2 dimensions.')
end

% shortcut
param = hamilt.pot;

% Atomic distances
R{1} = space.dvr.grid_ND{1};
R{2} = space.dvr.grid_ND{2};
R{3} = sqrt(R{1}.^2 + R{2}.^2 - 2*R{1}.*R{2}*cos(hamilt.pot.theta));

%% Now calculate the components, i.e., the coulomb and exchange integrals
for n = 1:3
    Q{n} = param.d_e(n) / (4 + 4*param.s(n))  * ( (3 + param.s(n)) ...
           * exp(-2*param.a(n)*(R{n} - param.r_e(n))) ...
           - (2 + 6*param.s(n)) * exp(-param.a(n)*(R{n} - param.r_e(n))) );

    J{n} = param.d_e(n) / (4 + 4*param.s(n)) * ( (1 + 3*param.s(n)) ...
           * exp(-2*param.a(n)*(R{n} - param.r_e(n))) ...
           - (6 + 2*param.s(n)) * exp(-param.a(n)*(R{n} - param.r_e(n))) );
end

%% And put them together for the complete LEPS surface
hamilt.pot.grid_ND{1,1} = Q{1} + Q{2} + Q{3} - ...
        sqrt( 0.5 * ((J{1} - J{2}).^2 + (J{1} - J{3}).^2 + (J{2} - J{3}).^2) );
