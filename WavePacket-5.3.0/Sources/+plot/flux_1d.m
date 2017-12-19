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
% see, e.g., Cohen-Tannoudji, Chapter III, Complement D, Eq. D-17 (p.239)
%
% normalization with 1/mass ?!?
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2013 Burkhard Schmidt
%
% see the README file for license details.

function flux_1d ( step, wide )
global expect hamilt info plots psi space

if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end

if hamilt.coupling.n_eqs~= 1
    util.error ('Cannot draw a flux density for >1 TDSEs')
end

%% Lower part: flux density
subplot ( 'Position', [1/w 1/h 7/w 3/h] );
hold off; plot( 1234567890, 1234567890 ); hold on;

% Calculate flux density
psi_mom = momentum(space.dof{1}, psi.dvr.grid_ND{1});
j = real ( conj(psi.dvr.grid_ND{1}) .* psi_mom );

% Curve plot of flux density
plot ( space.dvr.grid_1D{1}, ...
    j, ...
    'LineStyle', '-', ...
    'Color', plots.style.colors(1,:), ...
    'LineWidth', plots.style.line.thick )


%% Axes, labels, etc
axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max -1 1 ] )

set ( gca, 'XAxisLocation', 'bottom', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy)
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel('j(R)' )
    


%% Upper part: probability density
subplot ( 'Position', [1/w 5/h 7/w 3/h] );
hold off; plot( 1234567890, 1234567890 ); hold on;

% Get density from wavefunction
rho = abs   ( psi.dvr.grid_ND{1} ) .^2;

% Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
phi = angle ( psi.dvr.grid_ND{1} ) / (2*(pi+0.001)) + 1/2;

% Plot density and phase (with horizontal offset and base line)
if plots.density.energy.on
    offset = expect.ind.pot{1}(step);
else
    offset = 0;
end

plot.color ( space.dvr.grid_1D{1}, ...
    rho*plots.density.pot.delta/plots.density.rho_max.dvr, ...
    phi, ...
    plots.style.colors(1,:), ...
    plots.style.line.extrathick, ...
    offset, ...
    0 )

if plots.density.energy.on
    line ( [space.dof{1}.dvr_min space.dof{1}.dvr_max], ...
        [offset offset], ...
        'LineStyle', '-', ...
        'Color',     plots.style.colors(1,:), ...
        'LineWidth', plots.style.line.thin)
end

% Plot potential energy curve
if plots.density.energy.on
    plot ( space.dvr.grid_1D{1}, ...
        hamilt.pot.grid_ND{1,1}, ...
        'LineStyle', '-', ...
        'Color', plots.style.colors(1,:), ...
        'LineWidth', plots.style.line.thin )
end


% Axes, labels, etc
if plots.density.range.on == false
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
            plots.density.pot.min-plots.density.pot.delta/10 plots.density.pot.max+plots.density.pot.delta/10 ] )
else
    axis ( [ plots.density.range.x_min plots.density.range.x_max ...
            plots.density.pot.min-plots.density.pot.delta/10 plots.density.pot.max+plots.density.pot.delta/10 ] )
end

set ( gca, 'XAxisLocation', 'bottom', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy)
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
title ( {info.header1;info.header2} )
if plots.density.energy.on
    ylabel ( 'V(R)' )
else
    ylabel ( '\rho(R)' )
end





