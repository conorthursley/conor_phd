%------------------------------------------------------------------------------
%
% NOCl is described using two electronic states. The ground state potential
% is given as a superposition of harmonic oscillators as functions of the
% nuclear distances (which have to be calculated using the Jacobi coordinates...).
% The excited state is described by a polynomial expansion in some
% exponentiated Jacobi coordinates. The latter one needs quite a few
% terms.
%
% Parameters are the indices of the three dimensions in the setup
% - r_dof gives the NO bond distance (omitted for two dimensions)
% - R_dof gives the distance between Cl and NO CMS
% - c dof gives the bending angle 
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

function nocl

global hamilt space

%% Informational output

util.disp (' ')
util.disp ('*******************************************************')
util.disp ('Potential energy: NOCl                                 ')
util.disp ('Surface data taken from J. Chem. Phys 97:3199          ')
util.disp ('*******************************************************')


%% Check validity
if hamilt.coupling.n_eqs ~= 2
    util.error ('This potential only for two Schr�dinger equations')
end

if space.size.n_dim ~= 2 && space.size.n_dim ~= 3
    util.error ('This potential is only for two or three degrees of freedom')
end


%% Define fixed parameters

% The prefactors of the polynomial expansion of the excited state.
% See the reference for the meaning of these.
c_param = [ 0.0384816   0.0247875    0.0270933   0.00126791 0.00541285  0.0313629  0.0172449;
            0.00834237  0.00398713   0.00783319  0.0294887 -0.0154387  -0.0621984 -0.0337951;
            0.00161625 -0.00015633  -0.0189982  -0.00753297 0.00383665 -0.00758225 0.00904493;
           -0.0010101   0.000619148 -0.0149812  -0.0199722  0.00873013  0.0376118  0.0221523;
           -0.0003689   0.000164037 -0.00331809 -0.00567787 0.00268662  0.0134483  0.0084585;
           -0.0558666  -0.0276576    0.0934932  -0.0295638 -0.15436     0.0796119  0.135121;
            0.0582169   0.0384404    0.078114    0.185556  -0.0641656  -0.175976  -0.0104994;
            0.052285    0.0472724   -0.216008   -0.147775   0.349283    0.28458    0.00384449;
            0.0212609   0.0290597   -0.109124    0.0310445  0.262513   -0.250653  -0.369466;
            0.00334178  0.0039061   -0.0110452   0.0582029  0.0679524  -0.16459   -0.165337;
           -0.163186   -0.180535     0.04692     0.471673   0.403267   -0.718071  -0.761199;
           -0.0290674  -0.0136172   -0.108952   -1.68269   -1.2673      3.17648    2.92793;
            0.121228    0.202308     0.483613    1.29095   -0.174483   -2.4605    -1.36597;
            0.107233    0.115213    -0.366102    0.812662   1.76038    -1.19665   -1.77392;
            0.0232767   0.0304932   -0.19455    -0.0307517  0.539365    0.120203  -0.251289;
            0.0838975   0.198853    -0.0994766  -0.822409  -0.586006    1.17402    1.17378;
           -0.182047   -0.245637     0.130396    2.85439    2.44277    -5.36406   -5.22806;
           -0.227493   -0.470604    -0.670555   -1.66997    0.268677    3.71822    2.10678;
           -0.13635    -0.193843     0.626076   -1.55192   -3.22512     3.03851    4.01364;
           -0.0262554  -0.0391291    0.312858   -0.122063  -1.03112     0.28978    0.878604];

% Additional parameters for the excited state
a_param = [ 0.6816 -0.9123 0.4115 ];

% Scaling parameters for the exponentiated coordinates
alpha_param = 1.5;
beta_param = 1.1;

% Equilibrium positions of the Jacobi coordinates.
jacobi_re = [ 4.315 2.136 127.4 ];

% Masses needed for the conversion Jacobi => Bond length coordinates
mass_O = 17.9991610;
mass_N = 14.0030740048;

% oscillator strengths and anharmonic coupling for the ground state
k = [0.8987 0.0874 0.1137];
k23 = 0.0122;

% equilibrium bond distances
bl_re = [2.155 3.729 4.989];



%% Define the coordinates

% First, what are our Jacobi coordinates? For 2 dimensions, keep the NO
% distance fixed.
if space.size.n_dim == 2
    % we do not include the bond length
    Rgrid = space.dvr.grid_ND{hamilt.pot.R_dof};
    rgrid = ones(size(space.dvr.grid_ND{1})) * jacobi_re(2);
    cgrid = space.dvr.grid_ND{hamilt.pot.c_dof};
else
    Rgrid = space.dvr.grid_ND{hamilt.pot.R_dof};
    rgrid = space.dvr.grid_ND{hamilt.pot.r_dof};
    cgrid = space.dvr.grid_ND{hamilt.pot.c_dof};
end

% Exponentiated forms of the Jacobi coordinates for the excited state
qd = 1 - exp( -alpha_param * ( Rgrid - jacobi_re(1)) );
qv = rgrid - jacobi_re(2);
qt = exp(-beta_param * cgrid) - exp(-beta_param * cosd(jacobi_re(3)) );

% Bond length coordinates needed for the ground state
bl_r{1} = rgrid;
bl_r{2} = sqrt( (mass_O*rgrid / (mass_O+mass_N)).^2 + Rgrid.^2  ...
          + 2 * Rgrid .* rgrid .* cgrid * (mass_O / (mass_O+mass_N)) );
bl_r{3} = sqrt( (mass_N*rgrid/ (mass_O+mass_N) ).^2 + Rgrid.^2  ...
          - 2 * Rgrid .* rgrid .* cgrid * (mass_N / (mass_O+mass_N)) );



%% Create the potentials

% Ground state potential
hamilt.pot.grid_ND{1,1} = k23 * (bl_r{2}-bl_re(2)) .* (bl_r{3}-bl_re(3));
for d = 1:3
    hamilt.pot.grid_ND{1,1} = hamilt.pot.grid_ND{1,1} ...
                              + 0.5*k(d)*(bl_r{d}-bl_re(d)).^2;
end


% Excited state potential
hamilt.pot.grid_ND{2,2} = zeros(size(space.dvr.grid_ND{1}));

for l1 = 0:3
    for l2 = 0:4
        for l3 = 0:6
            hamilt.pot.grid_ND{2,2} = hamilt.pot.grid_ND{2,2} ...
                    + c_param(5*l1 + l2 + 1, l3 + 1) * qv.^l1 .* qd.^l2 .* qt.^l3;
        end
    end
end

hamilt.pot.grid_ND{2,2} = hamilt.pot.grid_ND{2,2} .* (1 - qd);

for l1 = 2:4
    hamilt.pot.grid_ND{2,2} = hamilt.pot.grid_ND{2,2} + a_param(l1-1) .* qv.^l1;
end


