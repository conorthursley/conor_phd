% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Ulf Lorenz
%
% see the README file for license details.

function H2O
global hamilt space

util.disp (' ')
util.disp ('****************************************************')
util.disp ('Semi-empirical calculation for the PES of water.    ')
util.disp ('The bound state is modeled after some Morse         ')
util.disp ('oscillators, the excited state based on ab-initio   ')
util.disp ('calculations. Fixed bending angle of 104.52 degrees.')
util.disp ('See J. Chem. Phys. 97:8285,( J. Chem. Phys. 88:129 )')
util.disp ('                                                    ')
util.disp ('Ground and excited state are joined such that       ')
util.disp ('dR(R1 = R_OH, R2 -> infinity) = 0                   ')
util.disp ('****************************************************')

%% Check validity
if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential only for two Schrödinger equations')
end

if space.size.n_dim ~= 2
    util.error ('This potential is for 2 degrees of freedom only')
end

%% Some variables we use throughout the code

gamma = 104.52; % constant HOH bending angle

R1 = space.dvr.grid_ND{1};      % O-H distance
R3 = space.dvr.grid_ND{2};      % other O-H distance
R2 = sqrt(R1.^2 + R3.^2 - 2*cosd(gamma) * R1.* R3);   % H-H distance

%% Ground state is just two Morse oscillators with a coupling term

% Note that we normalize the three separate atoms to E=0 in accordance
% with the paper by Engel and Schinke

D      = 0.2092;% dissociation energy  O-H
alpha  = 1.1327;% range parameter      O-H
r0     = 1.81;  % equilibrium distance O-H
F      = -6.76e-3;  % coupling coefficient
b      = 1;     % range parameter for coupling term

hamilt.pot.grid_ND{1,1} = D * (1 - exp(-alpha*(R1-r0))).^2 - D...
        + D * (1 - exp(-alpha*(R3-r0))).^2 - D...
        + F ./ (1 + exp(b .* ( (R1 - r0) + (R3 - r0) ))) .* (R1 - r0) .* (R3 - r0);


%% Excited state is a cumbersome fit to analytic data with asymptotic corrections

% Conversion becomes annoying after some time, so we convert the distances to
% Angstroem, and the final energy from eV to atomic units.
R1 = R1 * 0.5292;
R2 = R2 * 0.5292;
R3 = R3 * 0.5292;

S1 = R1 - 1;
S2 = R2 - 1;
S3 = R3 - 1;

% Parameters for the fit
delta = 0.25;   % range parameter for cutoff function
D_OH  = 4.621;  % dissociation energy
b_OH  = 2.294;  % range parameter
R_OH  = 0.971;  % equilibrium distance

a = [1.25 0.75 1.25];
c = [-0.9356263e0 -0.3399239e1 -0.3792337e1  0.3522993e1  0.2050857e1 ...
      0.3329077e0  0.6468300e0 -0.5144771e2 -0.4904665e1  0.2443589e2 ...
      0.7633727e0  0.1058277e2  0.7366955e1  0.6432985e2  0.8992322e1 ...
     -0.8969034e2  0.6452924e2 -0.5389425e2 -0.5909348e0 -0.1213870e2 ...
      0.6569770e1  0.5978727e2 -0.2984268e2  0.3034212e2 -0.2996398e2 ...
      0.2434441e2 -0.6282335e2  0.1819046e3  0.1035898e2 -0.3939324e2 ...
      0.3969865e1  0.1216507e2  0.5142234e1 -0.4749178e1 -0.1532321e1 ...
      0.9219124e1 -0.7433534e1  0.1573899e3  0.3329937e1 -0.3565261e1 ...
     -0.1541691e3 -0.1747232e1  0.5210177e2  0.1665088e3 -0.2398758e2 ...
     -0.6999091e2  0.1988446e2  0.4042948e2 -0.9719855e1  0.1853989e1];
% note that there is a typo in polynomial 10; the potential must
% be symmetric under exchange S1 <-> S3
power = { [0 0 0], [1 0 0], [0 1 0], [2 0 0], [0 2 0], ...
          [1 1 0], [1 0 1], [3 0 0], [0 3 0], [1 0 2], ...
          [1 2 0], [2 1 0], [1 1 1], [4 0 0], [0 4 0], ...
          [3 0 1], [2 0 2], [3 1 0], [2 2 0], [1 3 0], ...
          [1 2 1], [2 1 1], [5 0 0], [4 0 1], [3 0 2], ...
          [4 1 0], [3 1 1], [2 1 2], [3 2 0], [2 2 1], ...
          [2 3 0], [1 3 1], [1 4 0], [0 5 0], [6 0 0], ...
          [5 0 1], [4 0 2], [3 0 3], [5 1 0], [4 1 1], ...
          [3 1 2], [4 2 0], [3 2 1], [2 2 2], [3 3 0], ...
          [2 3 1], [2 4 0], [1 4 1], [1 5 0], [0 6 0]  };

% Now calculate the analytic fit.
% First the Morse terms, then the polynomial terms, and finally the tanh terms
% before we put it together.
pot1  = D_OH * (1 - exp(-b_OH * (R1 - R_OH))).^2 - D_OH;
pot3  = D_OH * (1 - exp(-b_OH * (R3 - R_OH))).^2 - D_OH;
morse = pot1 .* (1 - exp(-delta*R2.^2)) .* (1 - exp(-delta*R3.^2)) ...
      + pot3 .* (1 - exp(-delta*R2.^2)) .* (1 - exp(-delta*R1.^2));

poly = 0;
for ii = 1:50
    summand = c(ii) * S1.^power{ii}(1) .* S2.^power{ii}(2) .* S3.^power{ii}(3);
    if power{ii}(1) ~= power{ii}(3)
        summand = summand + c(ii) * S1.^power{ii}(3) .* S2.^power{ii}(2) .* S3.^power{ii}(1);
    end
    poly = poly + summand;
end

pot_tan = (1 - tanh(a(1) * S1)) .* (1 - tanh(a(2) * S2)) .* (1 - tanh(a(3)*S3));

pot_fit  = morse + poly.*pot_tan;

% It is not over yet. We have to construct an asymptotic potential for
% larger distances, which boils down to an analytic series again.

Rbar = 3.5*0.5292;  % new distance factor
A    = [-0.4579 2.6397 0.6379 0.134 2.7032 0.6428]; % constants

Rmax = R1;
Rmax(R3 > R1) = R3(R3>R1);

Rmin = R1;
Rmin(R3 < R1) = R3(R3<R1);

D_asy = D_OH + A(1) * exp( -A(2)*(Rmax - Rbar) - A(3)*(Rmax - Rbar).^2);
b_asy = b_OH + A(4) * exp( -A(5)*(Rmax - Rbar) - A(6)*(Rmax - Rbar).^2);

pot_asy = D_asy .* (1 - exp( -b_asy.*(Rmin-R_OH))).^2 - D_asy;

% Finally, join the fit and the asymptotic potential. First, we define the
% functions for joining, and perform the join thereafter

% omegas are defined in three regions. We first define them in the join region,
% and cut them off afterwards. At this point, of course, the paper uses atomic
% units again. For some reason, after I have seen their original Fortran code,
% I am not surprised.
omega_fit = 0.5 + 0.5 * cos(pi/2 * (Rmax - Rbar)/0.5292);
omega_asy = 0.5 - 0.5 * cos(pi/2 * (Rmax - Rbar)/0.5292);

omega_fit(Rmax < Rbar) = 1;
omega_asy(Rmax < Rbar) = 0;

omega_fit(Rmax > 5.5*0.5292) = 0;
omega_asy(Rmax > 5.5*0.5292) = 1;

% Now join the functions
hamilt.pot.grid_ND{2,2} = omega_fit .* pot_fit + omega_asy .* pot_asy;
hamilt.pot.grid_ND{2,2} = hamilt.pot.grid_ND{2,2} / 27.212; % conversion eV -> Hartree


%% Last thing to do: shift the excited state according to the recipe
%
% For large distances, the excited state uses some sort of Morse potential
% with different parameters. Right now, the two states are joined such that
% E_0(H + H + O) = E_1(H + H + O) = 0.
% However, this gives excitation energies that are significantly away from
% experimental data (53000 cm^-1 from the groundstate vs. >70000 cm^-1 at
% a single point for our data). To improve on this, we instead join the two 
% states such that
% E_0(OH + H) = E_1(OH + H), i.e. where one bond length is exactly at
% equilibrium distance, the other is infinity. It turns out that the 
% Shift is then just the difference of the dissociation energies.
hamilt.pot.grid_ND{2,2} = hamilt.pot.grid_ND{2,2} + (D_OH / 27.212) - D;
