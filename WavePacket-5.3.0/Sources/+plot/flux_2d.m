%--------------------------------------------------------------------------
%
% Calculate and visualize flux density (1D)
%
% J = -i*hbar/2 * [       psi^* * grad (psi) - psi * grad (psi^*) ]
% 
%   =    hbar * Im [      psi^* * grad (psi) ]
%
%   =    hbar * Re [ -i * psi^* * grad (psi) ]
%
%   =           Re [      psi^* *   P  (psi) ]
%
% see Cohen-Tannoudji, Chapter III, Complement D, Eq. D-17 (p.239)
%
% normalization with 1/mass ?!?
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.


function flux_2d ( step, wide )

global hamilt info plots psi space

if hamilt.coupling.n_eqs~= 1
    util.error ('Cannot draw a flux density for >1 TDSEs')
end

if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end
subplot ( 'Position', [1/w 1/h 7/w 7/h] );
hold off; plot( 1234567890, 1234567890 ); hold on;

% Plot dotted contours of total energy functions
if plots.density.energy.on
    contour ( space.dvr.grid_ND{1} , ...
        space.dvr.grid_ND{2}, ...
        hamilt.pot.grid_ND{1,1}, ...
        linspace( plots.density.pot.min, plots.density.pot.max, plots.density.contour.nlev(2) ), ...
        'LineStyle', ':', ...
        'LineWidth', plots.style.line.thin, ...
        'Color',     plots.style.colors(2,:)  );
end


lower = plots.density.rho_max.dvr / plots.density.contour.nlev(1);
upper = plots.density.rho_max.dvr;
if ~isempty(plots.density.contour.min)
    lower = plots.density.contour.min;
end
if ~isempty(plots.density.contour.max)
    upper = plots.density.contour.max;
end

% Contour plot of wavefunction
contour ( space.dvr.grid_ND{1},    ...
    space.dvr.grid_ND{2}, ...
    abs(psi.dvr.grid_ND{1}).^2, ...
    linspace(lower, upper,  plots.density.contour.nlev(1)), ...
    'LineStyle', '--', ...
    'LineWidth', plots.style.line.thin, ...
    'Color',     plots.style.colors(3,:)   );

% Quiver plot of flux density
psi_mom = momentum(space.dof{1}, psi.dvr.grid_ND{1});
j1 = real ( conj(psi.dvr.grid_ND{1}) .* psi_mom );
psi_mom = momentum(space.dof{2}, psi.dvr.grid_ND{1});
j2 = real ( conj(psi.dvr.grid_ND{1}) .* psi_mom );
quiver ( ...
    space.dvr.grid_ND{1}, ...
    space.dvr.grid_ND{2}, ...
    j1, ...
    j2, ...
    'LineStyle', '-', ...
    'LineWidth', plots.style.line.thick, ...
    'Color',     plots.style.colors(1,:)  );

if plots.density.range.on == false
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
            space.dof{2}.dvr_min space.dof{2}.dvr_max ] )
else
    range = plots.density.range;
    axis( [ range.x_min range.x_max range.y_min range.y_max ] )
end
set ( gca, 'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
title ( {info.header1;info.header2} )
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )

