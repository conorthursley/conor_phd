%--------------------------------------------------------------
%
% Visualize wavepacket in 2 dimensions 
% either in position (DVR) or in momentum (FBR) space 
% using surface plots
%
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function surface_2d ( step, wide )
global plots

% Create subplot of the right size
if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end
subplot ( 'Position', [1/w 1/h 7/w 7/h] );

% Either position (DVR) or in momentum (FBR) space 
if strcmpi ( plots.density.representation,'dvr' )
    density_dvr ( step )
elseif strcmpi ( plots.density.representation,'fbr' )
    density_fbr ( step )
end

%------------------------------------------------------------
% Plot 2-dim DVR densities and potential energy surfaces
%------------------------------------------------------------
function density_dvr ( step )
global expect hamilt info plots psi space

% Loop over wavefunctions/densities
for m=1:hamilt.coupling.n_eqs
    
    % Surface plots of potentials with densities color coded
    if plots.density.energy.on
        if expect.ind.pop{m}(step)>expect.min_pop
            surf ( space.dvr.grid_ND{1}, ...
                space.dvr.grid_ND{2}, ...
                hamilt.pot.grid_ND{m,m}, ...
                abs(psi.dvr.grid_ND{m}).^2 );
        else
            surf ( space.dvr.grid_ND{1}, ...
                space.dvr.grid_ND{2}, ...
                hamilt.pot.grid_ND{m,m} );
            
        end
        
        % Surface plots of densities only (and color-coded phases)
    else
        
        if expect.ind.pop{m}(step)>expect.min_pop
            
            % Normalized densities: Set maximum to + 1
            % Shift individual densities vertically
            rho = abs   ( psi.dvr.grid_ND{m} ) .^2;
            rho = min ( rho/plots.density.rho_max.dvr, 1 );
            rho = rho + 2*m - 2;

            % Dirty trick: Extend z-range from to -1 to 2*n-1
            rho( 1,  1) = 2*hamilt.coupling.n_eqs-1;
            rho(end,end) = -1;

            % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
            phi = angle ( psi.dvr.grid_ND{m} ) / (2*(pi+0.001)) + 1/2;
            phi(rho<0.01) = 0;
            
            % Surface plot of densities, with phases as colors
            surf ( space.dvr.grid_ND{1}, space.dvr.grid_ND{2}, rho, phi );

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
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )

if plots.density.range.on == false
    xyrange = [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
             space.dof{2}.dvr_min space.dof{2}.dvr_max ];
else
    r = plots.density.range;
    xyrange = [r.x_min r.x_max r.y_min r.y_max];
end

if plots.density.surface.color.on
    caxis( [plots.density.surface.color.min plots.density.surface.color.max] );
end

if plots.density.energy.on
    zrange = [ plots.density.pot.min-plots.density.pot.delta/10 plots.density.pot.max+plots.density.pot.delta/10 ];
    axis ( cat(2, xyrange, zrange) );
    zlabel ( 'Potential energy' )
    if hamilt.coupling.n_eqs==1
        zlabel ( 'V(R)' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi' )
            zlabel ( 'V_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia' )
            zlabel ( 'V_{dia}(R)' )
        end
    end
else
    zrange = [ -1 -1+2*hamilt.coupling.n_eqs ];
    axis ( cat(2, xyrange, zrange) );
    if hamilt.coupling.n_eqs==1
        zlabel ( '\rho(R) ' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi' )
            zlabel ( '\rho_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia' )
            zlabel ( '\rho_{dia}(R)' )
        end
    end
end

%------------------------------------------------------------
% Plot 2-dim FBR densities and kinetic energy surfaces
%------------------------------------------------------------
function density_fbr ( step )       
global expect hamilt info plots psi space


% Loop over wavefunctions/densities
for m=1:hamilt.coupling.n_eqs

    fbr = dvr2fbr(space.dof{1}, psi.dvr.grid_ND{m});
    fbr = dvr2fbr(space.dof{2}, fbr);

    % Surface plots of potentials with color coded densities 
    if plots.density.energy.on && isa(space.dof{1}, 'grid.fft') ...
        && isa(space.dof{2}, 'grid.fft')
        if expect.ind.pop{m}(step)>expect.min_pop
            surf ( space.fbr.grid_ND{1}, ...
                   space.fbr.grid_ND{2}, ...
                space.dof{1}.kin + space.dof{2}.kin, ...
                abs(fbr).^2 );
        else
            
            if m==1
                surf ( space.fbr.grid_ND{1}, ...
                       space.fbr.grid_ND{2}, ...
                       space.dof{1}.kin + space.dof{2}.kin );
            end
            
        end
        
    % Surface plots of densities only (and color-coded phases)
    else

        if expect.ind.pop{m}(step)>expect.min_pop

            % Normalized densities: Set maximum to + 1
            % Shift individual densities vertically
            rho = abs   ( fbr ) .^2;
            rho = min ( rho/plots.density.rho_max.fbr, 1 );
            rho = rho + 2*m - 2;

            % Dirty trick: Extend z-range from to -1 to 2*n-1
            rho( 1,  1 ) = 2*hamilt.coupling.n_eqs-1;
            rho(end,end) = -1;

            % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
            phi = angle ( fbr ) / (2*(pi+0.001)) + 1/2;
            phi(rho<0.01) = 0;
            
            % Surface plot of densities, with phases as colors
            surf ( space.fbr.grid_ND{1}, space.fbr.grid_ND{2}, rho, phi );

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
xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
ylabel ( [ 'P_{', space.dof{2}.label, '}' ] )

if plots.density.surface.color.on
    caxis( [plots.density.surface.color.min plots.density.surface.color.max] );
end

if plots.density.energy.on && isa(space.dof{1}, 'grid.fft') ...
        && isa(space.dof{2}, 'grid.fft')
    kin = space.dof{1}.kin + space.dof{2}.kin;
    maxkin = max(kin(:));
    axis ( [ space.dof{1}.fbr_min space.dof{1}.fbr_max ...
             space.dof{2}.fbr_min space.dof{2}.fbr_max ...
             0 maxkin ] )
    zlabel ( 'T(P)' )
else
    axis ( [ space.dof{1}.fbr_min space.dof{1}.fbr_max ...
             space.dof{2}.fbr_min space.dof{2}.fbr_max ...
             -1 -1+2*hamilt.coupling.n_eqs ] )
    if hamilt.coupling.n_eqs==1
        zlabel ( '\rho(P) ' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi' )
            zlabel ( '\rho_{adi}(P)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia' )
            zlabel ( '\rho_{dia}(P)' )
        end
    end
end


