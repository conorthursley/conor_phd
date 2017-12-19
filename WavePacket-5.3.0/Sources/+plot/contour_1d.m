%--------------------------------------------------------------------------
%
% Visualize Wigner transform of 1-dim wavepacket in phase space 
% using contour plots with marginals representing position
% and momentum densities. Alternatively, the density can also be
% drawn in DVR space.
%
% Compose animated figure from four subplots (see included subfunctions)
% containing densities (position, momentum, phase space) and respective
% energy functions (potential, kinetic, total)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2009 Burkhard Schmidt
%
% see the README file for license details.

function contour_1d ( step, wide )
global plots

if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end

% Without marginals
if ~plots.density.marginals.on

    % Wigner quasi-density, "mean trajectory" and total energy
    subplot ( 'Position', [1/w 1/h 7/w 7/h] )
    if ~plots.density.hold; hold off; end 
    plot( 1234567890, 1234567890 ); 
    hold on;
    density_wig ( step )

% With marginals
else
    
    % Upper left: Wigner quasi-density, "mean trajectory" and total energy
    subplot ( 'Position', [1/w 3/h 5/w 5/h] )
    if ~plots.density.hold; hold off; end 
    plot( 1234567890, 1234567890 );
    hold on;
    density_wig ( step )

    % Upper right: Momentum density and kinetic energy curve
    subplot ( 'Position', [6/w 3/h 2/w 5/h] )
    if ~plots.density.hold; hold off; end
    plot( 1234567890, 1234567890 ); 
    hold on;
    density_fbr  ( step )

    % Lower left: Position density and potential energy curve
    subplot ( 'Position', [1/w 1/h 5/w 2/h] )
    if ~plots.density.hold; hold off; end
    plot( 1234567890, 1234567890 );
    hold on;
    density_dvr ( step )

    % Lower right: Color mapping
    subplot ( 'Position', [6/w 1/h 2/w 2/h] )
    if ~plots.density.hold; hold off; end
    plot( 1234567890, 1234567890 );
    hold on;
    my_colors

end

%-----------------------------------------------------------
% Plot Wigner quasi-density and total energy surface
%-----------------------------------------------------------
function density_wig ( step )    
global expect hamilt info plots psi space uncert

% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs
    
    % Plot "mean trajectory" as solid curve with error bars in phase space: 
    % Expectation value of momentum versus expectation value of position
    if plots.density.expect.on
      plots.mask.ind{m} = find ( expect.ind.pop {m}(:) > expect.min_pop );
      errorbar ( expect.ind.dvr{m}(plots.mask.ind{m},1), ...
           expect.ind.fbr{m}(plots.mask.ind{m},1), ...
           uncert.ind.fbr{m}(plots.mask.ind{m},1), ...
          'LineStyle', '-', ...
          'LineWidth', plots.style.line.thick, ...
           'Color',     plots.style.colors(m,:) )
    end
       
    % Plot dotted contours of total enery functions
    if plots.density.energy.on
        contour ( space.dvr.grid_1D{1} , ...
                  -space.fbr.grid_1D{1}/2, ...
                  hamilt.tef.grid{m}(:,:), ...
                  linspace( plots.density.tef.min, plots.density.tef.max, plots.density.contour.nlev(2) ), ...
                  'LineStyle', ':', ...
                  'LineWidth', plots.style.line.thin, ...
                  'Color',     plots.style.colors(m,:)  );
    end
    
    % If population not too small, plot solid contours of Wigner functions
    if expect.ind.pop{m}(step)>expect.min_pop
        lower = -psi.wig.rho_max;
        upper = psi.wig.rho_max;
        if ~isempty(plots.density.contour.min)
            lower = plots.density.contour.min;
        end
        if ~isempty(plots.density.contour.max)
            upper = plots.density.contour.max;
        end

        contour ( space.dvr.grid_1D{1},    ...
                 -space.fbr.grid_1D{1}/2, ...
                 psi.wig.grid{m}(:,:), ...
                 linspace(lower, upper, plots.density.contour.nlev(1)), ...   % use even number of contours to avoid zero!
                 'LineStyle', '-', ...
                 'LineWidth', plots.style.line.thin, ...
                 'Color',     plots.style.colors(m,:)   );

    end

end

if plots.density.range.on == false
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
            space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2 ] )
else
    range = plots.density.range;
    axis ( [range.x_min  range.x_max  range.y_min  range.y_max ] )
end

set ( gca, 'XAxisLocation', 'top', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
xlabel ( {info.header1;info.header2} )
ylabel ( [ 'P_{', space.dof{1}.label, '}' ] )

% Negative imaginary potential (as absorbing boundary conditions)
if ~isempty(hamilt.nip.grid_ND)
    if hamilt.nip.min(1) > space.dof{1}.dvr_min
        line ( [ hamilt.nip.min(1) hamilt.nip.min(1)], ...
               [space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.max(1) < space.dof{1}.dvr_max
        line ( [ hamilt.nip.max(1) hamilt.nip.max(1)], ...
               [space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
end

%-------------------------------------------------
% Plot position density and potential energy curve
%-------------------------------------------------
function density_dvr ( step )
global expect hamilt plots psi space 


%% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs

    % If population not too small
    if ( expect.ind.pop{m}(step) > expect.min_pop )
        
        % Get density from wavefunction
        rho = abs   ( psi.dvr.grid_ND{m} ) .^2;

        % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
        phi = angle ( psi.dvr.grid_ND{m} ) / (2*(pi+0.001)) + 1/2;

        % Plot density and phase (with horizontal offset and base line)
        if plots.density.energy.on
            offset = expect.ind.pot{m}(step) + expect.ind.kin{m}(step);
        else
            offset = 0;
        end
        
        plot.color ( space.dvr.grid_1D{1}, ...
                     rho*plots.density.pot.delta/plots.density.rho_max.dvr, ...
                     phi, ...
                     plots.style.colors(m,:), ...
                     plots.style.line.thick, ...
                     offset, ...
                     0 )
       if plots.density.energy.on
            line ( [space.dof{1}.dvr_min space.dof{1}.dvr_max], ...
                   [offset offset], ...
                   'LineStyle', '-', ...
                   'Color',     plots.style.colors(m,:), ...
                   'LineWidth', plots.style.line.thin)
        end
        
    end

    % Plot potential energy curve
    if plots.density.energy.on
        plot ( space.dvr.grid_1D{1}, ...
               hamilt.pot.grid_ND{m,m}, ...
               'LineStyle', '-', ...
               'Color', plots.style.colors(m,:), ...
               'LineWidth', plots.style.line.thin )
    end
    
end

%% Axes, labels, etc
if plots.density.range.on == false
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
            plots.density.pot.min-plots.density.pot.delta/10 plots.density.pot.max+plots.density.pot.delta/10 ] )
else
    axis ( [ plots.density.range.x_min plots.density.range.x_max ...
            plots.density.pot.min-plots.density.pot.delta/10 plots.density.pot.max+plots.density.pot.delta/10 ] )
end

set ( gca, 'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy)
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )

if plots.density.energy.on

    if hamilt.coupling.n_eqs==1
        ylabel ( 'V(R)' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi')
            ylabel ( 'V_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia')
            ylabel ( 'V_{dia}(R)' )
        end
    end
    
else
    
    if hamilt.coupling.n_eqs==1
        ylabel ( '\rho(R)' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi')
            ylabel ( '\rho_{adi}(R)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia')
            ylabel ( '\rho_{dia}(R)' )
        end
    end

end

% Negative imaginary potential (as absorbing boundary condition)
if ~isempty(hamilt.nip.grid_ND)
    if hamilt.nip.min(1) > space.dof{1}.dvr_min
        line ( [hamilt.nip.min(1) hamilt.nip.min(1)], ...
               [plots.density.pot.min-plots.density.pot.delta/10 plots.density.pot.max+plots.density.pot.delta/10], ...
               'LineStyle', '--', ...
               'Color',     'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.max(1) < space.dof{1}.dvr_max
        line ( [hamilt.nip.max(1) hamilt.nip.max(1)], ...
               [plots.density.pot.min-plots.density.pot.delta/10 plots.density.pot.max+plots.density.pot.delta/10], ...
               'LineStyle', '--', ...
               'Color',     'k', ...
               'LineWidth', plots.style.line.thin)
    end
end

%----------------------------------------------------------
% Plot momentum density and kinetic energy curve (inverted)
%----------------------------------------------------------
function density_fbr ( step ) 
global expect hamilt plots psi space 

%% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs

    % If population not too small
    if ( expect.ind.pop{m}(step) > expect.min_pop )

        psi.fbr.grid_ND{m} = dvr2fbr(space.dof{1}, psi.dvr.grid_ND{m});
        
        % Get density from wavefunction
        rho = abs ( psi.fbr.grid_ND{m}(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4) ) .^2;

        % Get phase of wavefunction: Map interval [-pi,pi] into [0,1]
        phi = angle ( psi.fbr.grid_ND{m}(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4) ) / (2*(pi+0.001)) + 1/2;

        % Plot density and phase (with horizontal offset and base line)
        if plots.density.energy.on
            offset = expect.ind.pot{m}(step) + expect.ind.kin{m}(step);
        else
            offset = 0;
        end
            
        plot.color ( space.fbr.grid_1D{1}(space.dof{1}.n_pts/4+1 : 3*space.dof{1}.n_pts/4), ...
                      rho*plots.density.kin.delta/plots.density.rho_max.fbr, ...
                      phi, ...
                      plots.style.colors(m,:), ...
                      plots.style.line.thick, ...
                      offset, ...
                      1 )
                  
        if plots.density.energy.on
               line ( [ offset        offset       ], ...
                      [space.dof{1}.fbr_min space.dof{1}.fbr_max], ...
                      'LineStyle', '-', ...
                      'Color',     plots.style.colors(m,:), ...
                      'LineWidth', plots.style.line.thin)
        end
                      
    end

end

%% Plot kinetic energy curve
if plots.density.energy.on
    plot ( space.dof{1}.kin, ...
           -space.fbr.grid_1D{1}, ...
           'Color',    'k', ...
           'LineWidth', plots.style.line.thin )
end
    
%% Axes, labels, etc
if plots.density.range.on == false
    axis ( [ plots.density.kin.min-plots.density.kin.delta/10 plots.density.kin.max+plots.density.kin.delta/10 ...
            space.dof{1}.fbr_min/2 space.dof{1}.fbr_max/2 ] )
else
    axis ( [ plots.density.kin.min-plots.density.kin.delta/10 plots.density.kin.max+plots.density.kin.delta/10 ...
            plots.density.range.y_min plots.density.range.y_max ] )
end

set ( gca, 'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
if plots.density.energy.on

        xlabel ( 'T(P)' )
    
else
    
    if hamilt.coupling.n_eqs==1
        xlabel ( '\rho(P)' )
    else
        if strcmpi ( hamilt.coupling.representation,'adi' )
            xlabel ( '\rho_{adi}(P)' )
        elseif strcmpi ( hamilt.coupling.representation,'dia' )
            xlabel ( '\rho_{dia}(P)' )
        end
    end

end
ylabel ( [ 'P_{', space.dof{1}.label, '}' ] )

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
           'YAxisLocation', 'right', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
xlabel ('phase [\pi]')
ylabel ('color map')

