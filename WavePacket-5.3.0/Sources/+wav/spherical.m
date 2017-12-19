%------------------------------------------------------------------------------
%
% This function creates associated Legendre polynomials in cos Theta
%
% Apart from normalization and missing azimuthal functions,
% these functions are identical to spherical harmonics.
%
% Optionally, dividing wavefunction by sqrt(sin(Theta))
% see Eq. (4) in doi:10.1103/PhysRev.A.91.022111
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2017 Burkhard Schmidt
%
% see the README file for license details.

function init_grid = spherical(dir)
global space psi

%% Set default values
if ~isfield(psi.dof{dir}, 'l')
        psi.dof{dir}.l = 0;
end
l = psi.dof{dir}.l;

if ~isfield(psi.dof{dir},'m')
    psi.dof{dir}.m = 0;
end
m = psi.dof{dir}.m;

if ~isfield(psi.dof{dir},'sqst')
    psi.dof{dir}.sqst = false;
end

% Check input
if l<0
    util.error ('quantum number l must not be negative')
end
if m<0
    util.error ('quantum number m must not be negative')
end
if m>l
    util.error ('quantum number m must not exceed l')
end


%% Some output
util.disp (' ')
util.disp ('*******************************************************')
util.disp ( ['Initial wavefunction for DOF :' int2str(dir)] )
util.disp ('   ' )
util.disp ('Spherical harmonic')
util.disp ( ['Angular momentum quantum number l : ' int2str(l)] )
util.disp ( ['Azimuthal quantum number        m : ' int2str(m)] )
util.disp ( ['Divide by sqrt(sin(Theta))        : ' int2str(psi.dof{dir}.sqst)] )


%% Set up the grid, using Matlab's built-in Legendre polynomials
X = cos ( space.dvr.grid_ND{dir} );
Y = legendre (l,X);
init_grid = Y (m+1,:)';

if psi.dof{dir}.sqst
    init_grid = init_grid .* sqrt( sin ( space.dvr.grid_ND{dir} ) );
end
