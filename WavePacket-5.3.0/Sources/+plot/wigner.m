%------------------------------------------------------------------------------
%
% Prepare Wigner plots of 1-dim wavepacket dynamics:
% Grid representations of energy and Wigner functions
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%
% see the README file for license details.

function wigner ( step )
global expect hamilt psi plots space

% Wigner plot only for 1D problems
if ~plots.density.on 
    return
end
if  space.size.n_dim ~= 1 ...
        || (~strcmp(plots.density.type, 'contour') && ~strcmp(plots.density.type, 'surface'))
    return
end

% Wigner transform for FFT grids only
if ~isa (space.dof{1}, 'grid.fft')
    util.disp ('Cannot make a Wigner plot of a 1D non-Fourier grid. Using "curve" instead')
    plots.density.type = 'curve';
    return
end

% Grid representation of total energy function(s)
if plots.density.energy.on && step==1

    % Make grid representation of total energy function compatible
    % with Wigner representation of wavefunction:
    % Use only inner part of kinetic energy grid and double size
    kin2 = space.dof{1}.kin;
    kin2 = kin2(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4);
    kin2 = [kin2 kin2]';
    kin2 = reshape ( kin2, space.dof{1}.n_pts, 1 );

    for m = 1:hamilt.coupling.n_eqs
        hamilt.tef.grid {m}  = ...
            repmat ( hamilt.pot.grid_ND{m,m}, 1, space.dof{1}.n_pts)' + ...
            repmat (        kin2            , 1, space.dof{1}.n_pts);

        hamilt.tef.grid{m}( hamilt.tef.grid{m}>hamilt.truncate.max ) = hamilt.truncate.max;
        hamilt.tef.grid{m}( hamilt.tef.grid{m}<hamilt.truncate.min ) = hamilt.truncate.min;

    end
end


% Perform Wigner transform of all wavefunctions
for m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step) > expect.min_pop
        psi.wig.grid{m} = ket.wigner ( psi.dvr.grid_ND{m}' );
    end
end
