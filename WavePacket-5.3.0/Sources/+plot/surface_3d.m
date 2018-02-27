%--------------------------------------------------------------
%
% Visualize wavepacket in 3 dimensions 
% so far: only in position (DVR)space 
% using Matlab's builtin isosurface plots
%
%--------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.

function surface_3d ( step, wide )
global hamilt info plots psi space
persistent x y z rhomax

if hamilt.coupling.n_eqs~= 1
    util.error ('Can not yet draw a 3D isosurface for coupled TDSEs')
end

% Create subplot of the right size
if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end
subplot ( 'Position', [1/w 1/h 7/w 7/h] );

% Get (position or momentum) density
if strcmpi ( plots.density.representation,'dvr' )
    rho = abs ( psi.dvr.grid_ND{1} ) .^2;
elseif strcmpi ( plots.density.representation,'fbr' )
    fbr = dvr2fbr(space.dof{1}, psi.dvr.grid_ND{1});
    fbr = dvr2fbr(space.dof{2}, fbr);
    fbr = dvr2fbr(space.dof{3}, fbr);
    rho = abs ( fbr ) .^2;
end

% Save x/y/z and find maximum of density (first step only)
if step==1
    % No need to arrange arrays x,y,z as if they were created with MESHGRID
    if strcmpi ( plots.density.representation,'dvr' )
        x = space.dvr.grid_ND{1};
        y = space.dvr.grid_ND{2};
        z = space.dvr.grid_ND{3};
    elseif strcmpi ( plots.density.representation,'fbr' )
        x = space.fbr.grid_ND{1};
        y = space.fbr.grid_ND{2};
        z = space.fbr.grid_ND{3};
    end
    rhomax = max(abs(rho(:)));
else
    cla %  deletes all visible graphics objects from the current axes
end

% Surface plot of densities
patch(isosurface ( x,y,z, rho, rhomax/10), ...
    'FaceColor', 'green', ...
    'EdgeColor', 'none');

% Specify view point in terms of azimuth and elevation
view (plots.density.surface.view(1),plots.density.surface.view(2))

% Lighting of surfaces
camlight;  
lighting PHONG;
if plots.density.surface.look(2)
    lightangle(plots.density.surface.lite(1),plots.density.surface.lite(2));
end

% Axes, labels, etc
if strcmpi ( plots.density.representation,'dvr' )
    xyrange = [ ...
        space.dof{1}.dvr_min space.dof{1}.dvr_max ...
        space.dof{2}.dvr_min space.dof{2}.dvr_max ...
        space.dof{3}.dvr_min space.dof{3}.dvr_max ];
elseif strcmpi ( plots.density.representation,'fbr' )
    xyrange = [ ...
        space.dof{1}.fbr_min space.dof{1}.fbr_max ...
        space.dof{2}.fbr_min space.dof{2}.fbr_max ...
        space.dof{3}.fbr_min space.dof{3}.fbr_max ];
end
axis (xyrange);

set ( gca, 'LineWidth',     plots.style.line.thick, ...
    'FontName',      plots.style.font.name,  ...
    'FontSize',      plots.style.font.large, ...
    'FontWeight',    plots.style.font.heavy )
title  ( {info.header1;info.header2} )
if strcmpi ( plots.density.representation,'dvr' )
    xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
    ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )
    zlabel ( [ 'R_{', space.dof{3}.label, '}' ] )
elseif strcmpi ( plots.density.representation,'fbr' )
    xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
    ylabel ( [ 'P_{', space.dof{2}.label, '}' ] )
    zlabel ( [ 'P_{', space.dof{3}.label, '}' ] )
end
