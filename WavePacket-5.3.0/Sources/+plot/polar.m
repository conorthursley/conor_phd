%--------------------------------------------------------------------------
%
% Polar plots of one-dimensional wavefunctions/densities on a periodic
% domain. The coordinate is the angle or its cosine, depending on the grid),
% the distance from the origin is the squared norm, the color is the phase.
% Since the plotting is a bit difficult to interpret, no potentials are
% ever plotted.
%
% Matlab offers a polar plot which we do not use for a couple of reasons.
% First, it does not really set up a polar plot, but just draws an "image"
% of a polar plot in an underlying x/y plot. You cannot really resize it etc.
% Second, it cannot be configured; zero angle is in the horizontal, which is
% very unintuitive at least for me. So it boils down to drawing a couple of
% lines that we can also do by hand.
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2009 Ulf Lorenz
%               2009 Burkhard Schmidt
%
% see the README file for license details.

function polar ( step, wide )
global space

if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end

if space.size.n_dim ~= 1
    util.error('Cannot draw a polar plot for >1 dimensions.')
end

% subplot ( 'Position', [1/w 3/h 7/w 5/h] );
subplot ( 'Position', [1/w 1/h 7/w 7/h] );
hold off; plot( 1234567890, 123456789); hold on;

density_polar ( step );

% subplot ( 'Position', [1/w 1/h 7/w 1/h] )
% my_colors

end

%-------------------------------------------------
% Plot position density
%-------------------------------------------------
function density_polar ( step )
global expect hamilt info plots psi space 


%% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs

    % If population not too small
    if ( expect.ind.pop{m}(step) > expect.min_pop )
        
        % Get density from wavefunction
        rho = abs   ( psi.dvr.grid_ND{m} ) .^2 / plots.density.rho_max.dvr;

        % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
        phi = angle ( psi.dvr.grid_ND{m} ) / (2*(pi+0.001)) + 1/2;

        % Get the angle. If we have a grid.legendre, we assume the
        % coordinate is the cosine of the angle, otherwise it is
        % the angle
        if isa(space.dof{1}, 'grid.legendre')
            theta = acos(space.dvr.grid_ND{1});
            % cover the entire range [0,2pi] instead of [0,pi]
            theta = cat(1, theta, 2*pi - flipud(theta));
            rho   = cat(1, rho, flipud(rho));
            phi   = cat(1, phi, flipud(phi));
        else
            theta = space.dvr.grid_ND{1};
            % HACK: add the first point of theta to the end, so that we
            % have a periodic line drawing
            theta = cat(1, theta, theta(1));
            rho   = cat(1, rho, rho(1));
            phi   = cat(1, phi, phi(1));
        end

        % convert position to Cartesian coordinates; note theta=0 is vertical
        x = rho .* sin(theta);
        y = rho .* cos(theta);

        % Draw colored curve; extra-thick!
        plot.color(x, y, phi, plots.style.colors(m,:), plots.style.line.extrathick, 0, 2)
    end
end

% Draw a few radial lines to guide the eye in steps of 30 degrees
for theta = linspace(0, pi, 7)
    line([-sin(theta) sin(theta)], [-cos(theta) cos(theta)], ...
            'LineStyle', ':', 'LineWidth', plots.style.line.thin);
end

% Draw a few circles to guide the eye in steps of 0.25 radius
for r = (0.25:0.25:1)
    angles = linspace(0,2*pi,100);
    plot(r*cos(angles), r*sin(angles), ':', 'LineWidth', plots.style.line.thin);
end

%% Axes, labels, etc
axis ( [ -1 1 -1 1] )

set ( gca, 'XAxisLocation','top', ...
           'XTickLabel',{}, ...
           'YTickLabel',{}, ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy)
xlabel ( {info.header1;info.header2} )
end


%----------------------------------------------------------
% Plot color mapping
%----------------------------------------------------------
function my_colors 
global plots

x = linspace ( -1, 1, 50 );
y = ones (size(x));
c = linspace ( 0, 1, 50 );
plot.color ( x, y, c, [1 0 0], plots.style.line.extrathick, 0, 0)
axis ( [-1 1 0.9 1.1] )
set ( gca, 'XTick',          -0.5:0.5:1, ...
           'YTick',         [-123456789 123456789], ... % Fake: suppress tick labels
           'YAxisLocation', 'left' , ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
xlabel ('phase [\pi]')
ylabel (['color';' map '])

end
