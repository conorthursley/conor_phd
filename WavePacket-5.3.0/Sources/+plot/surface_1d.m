%--------------------------------------------------------------
%
% Visualize Wigner transform of 1-dim wavepacket in phase space 
% using surface plots
%
% Animated graphics of (Wigner) quasi-density in phase space
% with marginal distributions (i.e. position/momentum densities)
%
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function surface_1d ( step, wide )
global expect hamilt info plots psi space 

if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end
subplot ( 'Position', [1/w 1/h 7/w 7/h] );

% Loop over wavefunctions/densities
for m=1:hamilt.coupling.n_eqs

    % Surface plot of energy functions with color coded Wigner functions
    if plots.density.energy.on
        if expect.ind.pop{m}(step)>expect.min_pop
            surf ( space.dvr.grid_1D{1}, ...
                space.fbr.grid_1D{1}/2, ...
                 hamilt.tef.grid{m}, ...
                 psi.wig.grid{m} );
        else
            surf ( space.dvr.grid_1D{1}, ...
                space.fbr.grid_1D{1}/2, ...
                hamilt.tef.grid{m} );
        end
        
    % Surface plots of Wigner functions
    else

        if expect.ind.pop{m}(step)>expect.min_pop

            % Normalize Wigner distribution: Set maximum to +/- 1
            % Shift individual quasi-densities vertically
            psi.wig.grid{m} = min( psi.wig.grid{m}/psi.wig.rho_max, 1 ) + 2*m - 2;

            if plots.density.marginals.on

                % Position densities in first row as marginal distribution
                rho_pos = abs ( psi.dvr.grid_ND{m} ).^2;
                rho_pos = min( rho_pos/plots.density.rho_max.dvr, 1 );
                psi.wig.grid {m}(1,:) = psi.wig.grid {m}(1,:) + rho_pos';

                % Take only inner part of momentum function and double the size
                psi.mom = ifft ( psi.dvr.grid_ND{m} ) * sqrt(2*pi) / ...
                            (space.fbr.grid_1D{1}(2) - space.fbr.grid_1D{1}(1));    % delta_p
                psi.mom = fftshift ( psi.mom );
                psi.mom = psi.mom ( space.dof{1}.n_pts/4+1:3*space.dof{1}.n_pts/4 );
                psi.mom = [psi.mom psi.mom]';
                psi.mom = reshape ( psi.mom, space.dof{1}.n_pts, 1 )/2;

                % Momentum densities in first column as marginal distribution
                rho_mom = abs ( psi.mom ).^2;
                rho_mom = min( rho_mom/plots.density.rho_max.fbr*4, 1 );
                psi.wig.grid {m}(:,1) = psi.wig.grid {m}(:,1) + rho_mom;

            end

            % Dirty trick: Extend z-range from to -1 to 2n-1
            psi.wig.grid {m}(1,1) = 2*hamilt.coupling.n_eqs-1;
            psi.wig.grid {m}(2,2) = -1;

            % Surface plot of phase space (Wigner) distribution
            surf ( space.dvr.grid_1D{1}, -space.fbr.grid_1D{1}/2, psi.wig.grid{m} );

        end

    end

    if m==1
        hold on
    end

end
hold off

% Specify view point in terms of azimuth and elevation
view (plots.density.surface.view(1),plots.density.surface.view(2))

% Select color map
colormap(plots.density.surface.color.map);

% Shading of surface(s)
if plots.density.surface.look(1)
    shading interp;
end

% Lighting of surfaces
if plots.density.surface.look(2)
    lightangle(plots.density.surface.lite(1),plots.density.surface.lite(2));
end

% Axes, labels, etc
set ( gca, 'LineWidth',     plots.style.line.thick, ...
    'FontName',      plots.style.font.name,  ...
    'FontSize',      plots.style.font.large, ...
    'FontWeight',    plots.style.font.heavy )
title  ( {info.header1;info.header2} )
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'P_{', space.dof{1}.label, '}' ] )

if plots.density.range.on == false
    xyrange = [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
            -space.dof{1}.fbr_max/2 space.dof{1}.fbr_max/2 ];
else
    r = plots.density.range;
    xyrange = [ r.x_min r.x_max r.y_min r.y_max ];
end

if plots.density.surface.color.on
    caxis( [plots.density.surface.color.min plots.density.surface.color.max] )
end

if plots.density.energy.on
    zrange = [ plots.density.tef.min-plots.density.tef.delta/10 plots.density.tef.max+plots.density.tef.delta/10 ];
    axis ( cat(2, xyrange, zrange) );
    if hamilt.coupling.n_eqs==1
        zlabel ( 'E(R,P)' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi')
            zlabel ( 'E_{adi}(R,P)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia')
            zlabel ( 'E_{dia}(R,P)' )
        end
    end
else
    zrange = [ -1 -1+2*hamilt.coupling.n_eqs ];
    axis ( cat(2, xyrange, zrange) );
    if hamilt.coupling.n_eqs==1
        zlabel ( '\rho(R,P)' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi')
            zlabel ( '\rho_{adi}(R,P)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia')
            zlabel ( '\rho_{dia}(R,P)' )
        end
    end
end


