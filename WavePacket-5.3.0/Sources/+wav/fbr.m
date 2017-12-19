%------------------------------------------------------------------------------
%
% This function creates an initial state as a "momentum", or rather FBR eigenstate,
% (single state or coherent superposition of states)
% where the exact nature of the state depends on the grid used in the DVR scheme.
% For FFT grids, it returns a plane wave, for Legendre polynomials, the result
% is a spherical harmonic, for Hermite grids it is a Harmonic oscillator eigenstate.
%
% However, since you have to specify that you want the n-th eigenfunction, this
% routine is somewhat awkward to use for plane waves, as their quantum numbers
% translate into momenta going from -pmax to +pmax. I.e. to get momentum 0, you
% have to specify something like "#grid points / 2".
%
% You can either set up a pure eigenstate by specifying the variable state or a
% superposition by specifying the array of coefficients in the variable coeffs.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2007-2008 Ulf Lorenz
%               2008 Burkhard Schmidt
%
% see the README file for license details.

function init_grid = fbr(dir)

global space psi

util.disp ( ' ' )
util.disp ( '*******************************************************' )
util.disp ( ['Initial wavefunction for DOF :' int2str(dir)]           )
util.disp ( 'grid eigenfunction.                                    ' )
if isfield(psi.dof{dir}, 'state')
    util.disp ( ['We pick the ' int2str(psi.dof{dir}.state) '-th eigenstate'])

    init_grid = zeros(size(space.fbr.grid_ND{dir}));
    init_grid( space.fbr.grid_ND{dir} == space.fbr.grid_1D{dir}(psi.dof{dir}.state)) = 1;

elseif isfield(psi.dof{dir}, 'coeffs')
    util.disp('Using user-specified coefficients')

    init_grid = zeros(size(space.fbr.grid_ND{dir}));
    for i = 1:length(psi.dof{dir}.coeffs)
        init_grid(space.fbr.grid_ND{dir} == space.fbr.grid_1D{dir}(i)) = psi.dof{dir}.coeffs(i);
    end

else
    util.err('Cannot set up initial wave function; parameters not supplied')
end

% Transform the result back to DVR space
init_grid = fbr2dvr(space.dof{dir}, init_grid);
