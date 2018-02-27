% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2008 Ulf Lorenz
%
% see the README file for license details.

function H3plus
global hamilt space

util.disp (' ')
util.disp ('***************************************************')
util.disp ('High quality MRCI calculation for potential energy ')
util.disp ('surface of  H3+                                    ')
util.disp ('See J. Chem. Phys. 84:891                          ')
util.disp ('                                                   ')
util.disp ('***************************************************')
util.disp ( [ 'angular coordinate : ' num2str(hamilt.pot.c_dof) ] )

%% Check validity
if hamilt.coupling.n_eqs ~= 1
    util.error ('This potential only for single Schrödinger equation')
end

if space.size.n_dim ~= 3
    util.error ('This potential is for 3 degrees of freedom only')
end

if ~isa (space.dof{hamilt.pot.c_dof}, 'grid.legendre')
    util.error ('The potential requires Jacobi coordinates.')
end


%% Transform variables

% The transformation goes in three steps.
% 1. We transform from Jacobi to Bond length coordinates
% 2. The bond length coordinates are transformed in an exponential form
% 3. Using these exponentiated coordinates, one can define the S3
%    symmetry-adapted deformation coordinates
%
% The potential is then given as a power series in the deformation
% coordinates.


% First, define the angle and the Jacobi lengths
if hamilt.pot.c_dof == 1
    ang = space.dvr.grid_ND{1};
    j1 = space.dvr.grid_ND{2};
    j2 = space.dvr.grid_ND{3};
elseif hamilt.pot.c_dof == 2
    ang = space.dvr.grid_ND{2};
    j1 = space.dvr.grid_ND{1};
    j2 = space.dvr.grid_ND{3};
else
    ang = space.dvr.grid_ND{3};
    j1 = space.dvr.grid_ND{1};
    j2 = space.dvr.grid_ND{2};
end

% Now go over to bond length coordinates
r12 = j1;
r13 = sqrt( j1.^2/4 + j2.^2 + 2 * j1/2 .* j2 .* ang );
r23 = sqrt( j1.^2/4 + j2.^2 - 2 * j1/2 .* j2 .* ang );

% Transform to the exponential form
beta = 1.3;
R_e  = 1.6501;
q12 = 1/beta * ( 1 - exp(-beta * (r12/R_e-1)) );
q13 = 1/beta * ( 1 - exp(-beta * (r13/R_e-1)) );
q23 = 1/beta * ( 1 - exp(-beta * (r23/R_e-1)) );

% Finally, transform to deformation coordinates
sa = (q12 + q13 + q23) / sqrt(3);
sx = (2*q12 - q23 - q13) / sqrt(6);
sy = (q23 - q13) / sqrt(2);

% The coordinates x/y can be expressed by magnitude time cos/sin angle
se = sqrt(sx.^2 + sy.^2);
phi = acos(sx ./ se);


%% Create the potential

% The potential is given as a power series expansion
%                     n   2m+3k
% V    = sum    V    S   S       cos (3k*phi)
%  nmk      nmk  nmk  a   e
%
% with n + 2m + 3k <= N. we use the expansion for N = 7.

V = zeros(8, 4, 3);

% V_nm0
V(:,:, 1) = [ 0      266219  44851  15068;
              130   -241851 -11820  37493;
              204603 131115  120688 0;
             -49832 -50919   23840  0;
              25002  50424   0      0;
             -2115   887     0      0;
              4346   0       0      0;
             -277    0       0      0];

% V_nm1
V(:,:,2) = [-6490 -3185   7605 0;
             88648 73273  0    0;
            -28688 104361 0    0;
             57028 0      0    0;
             9333  0      0    0;
             0     0      0    0;
             0     0      0    0;
             0     0      0    0];

% V_nm2
V(:, 1, 3) = [-339; -3238; 0; 0; 0; 0; 0; 0];

V = V * 1e-6;

% And finally, create the potential
hamilt.pot.grid_ND{1,1} = zeros(size(space.dvr.grid_ND{1,1}));

for n = 0:7
    for m = 0:3
        for k = 0:2
            hamilt.pot.grid_ND{1,1} = hamilt.pot.grid_ND{1,1} + ...
                            V(n+1,m+1,k+1) * sa.^n .* se.^(2*m+3*k) .* cos(3*k*phi);
        end
    end
end
