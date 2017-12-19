%--------------------------------------------------------------------------
%
% Visualize 2-dim wavepacket in DVR or in FBR space 
% using contour plots
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%
% see the README file for license details.


function contour_2d ( step, wide )
global plots

if wide % Wide format: 16:9
    w=16; h=09;
else % Square format: 9:9
    w=09; h=09;
end

% Without marginals
if ~plots.density.marginals.on

    % 2-dim contour plot
    subplot ( 'Position', [1/w 1/h 7/w 7/h] )
    if ~plots.density.hold; hold off; end 
    plot( 1234567890, 1234567890 ); 
    hold on;
    if strcmpi ( plots.density.representation,'dvr' )
        density_dvr ( step )
    elseif strcmpi ( plots.density.representation,'fbr' )
        density_fbr ( step )
    end

% With marginals
else
    
    % Upper left: 2-dim contour plot
    subplot ( 'Position', [1/w 3/h 5/w 5/h] )
    if ~plots.density.hold; hold off; end 
    plot( 1234567890, 1234567890 ); 
    hold on;
    if strcmpi ( plots.density.representation,'dvr' )
        density_dvr ( step )
    elseif strcmpi ( plots.density.representation,'fbr' )
        density_fbr ( step )
    end

    % Upper right: Projection on second coordinate
    subplot ( 'Position', [6/w 3/h 2/w 5/h] )
    if ~plots.density.hold; hold off; end 
    plot( 1234567890, 1234567890 ); 
    hold on;
    if strcmpi ( plots.density.representation,'dvr' )
        density_dvr2 ( step )
    elseif strcmpi ( plots.density.representation,'fbr' )
        density_fbr2 ( step )
    end

    % Lower left: Projection on first coordinate
    subplot ( 'Position', [1/w 1/h 5/w 2/h] )
    if ~plots.density.hold; hold off; end  
    plot( 1234567890,1234567890 ); 
    hold on;
    if strcmpi ( plots.density.representation,'dvr' )
        density_dvr1 ( step )
    elseif strcmpi ( plots.density.representation,'fbr' )
        density_fbr1 ( step )
    end

end

%------------------------------------------------------------
% Plot 2-dim position densities and potential energy surfaces
%------------------------------------------------------------
function density_dvr ( step )       
global expect hamilt info plots psi space

% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs
    
    % Plot "mean trajectory" as solid curve in position space
    % Expectation value of R_1 vs. R_2
    if plots.density.expect.on
      plots.mask.ind{m} = find ( expect.ind.pop{m} > expect.min_pop );
      plot ( expect.ind.dvr{m}(plots.mask.ind{m},1), ...
           expect.ind.dvr{m}(plots.mask.ind{m},2), ...
          'LineStyle', '-', ...
          'LineWidth', plots.style.line.thick, ...
           'Color',     plots.style.colors(m,:) )
    end
       
    % Plot dotted contours of total energy functions
    if plots.density.energy.on
        contour ( space.dvr.grid_ND{1} , ...
                  space.dvr.grid_ND{2}, ...
                  hamilt.pot.grid_ND{m,m}, ...
                  linspace( plots.density.pot.min, plots.density.pot.max, plots.density.contour.nlev(2) ), ...
                  'LineStyle', ':', ...
                  'LineWidth', plots.style.line.thin, ...
                  'Color',     plots.style.colors(m,:)  );
    end
    
    % If populations not too small, plot solid contours of densities
    if expect.ind.pop{m}(step)>expect.min_pop
        lower = plots.density.rho_max.dvr / plots.density.contour.nlev(1);
        upper = plots.density.rho_max.dvr;
        if ~isempty(plots.density.contour.min)
            lower = plots.density.contour.max;
        end
        if ~isempty(plots.density.contour.min)
            upper = plots.density.contour.max;
        end
        
        contour ( space.dvr.grid_ND{1},    ...
                  space.dvr.grid_ND{2}, ...
                  abs(psi.dvr.grid_ND{m}).^2, ...
                  linspace(lower, upper,  plots.density.contour.nlev(1)), ...   
                 'LineStyle', '-', ...
                 'LineWidth', plots.style.line.thin, ...
                 'Color',     plots.style.colors(m,:)   );
    end
    
end

if plots.density.range.on == false
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ...
            space.dof{2}.dvr_min space.dof{2}.dvr_max ] )
else
    range = plots.density.range;
    axis( [ range.x_min range.x_max range.y_min range.y_max ] )
end
set ( gca, 'XAxisLocation', 'top', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
title ( {info.header1;info.header2} )
if ~plots.density.marginals.on
    xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
end
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )

% Absorbing boundary conditions
if ~isempty(hamilt.nip.grid_ND)
    if hamilt.nip.min(1) > space.dof{1}.dvr_min % Left
        line ( [ hamilt.nip.min(1) hamilt.nip.min(1) ], ...
               [ space.dof{2}.dvr_min space.dof{2}.dvr_max ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.max(1) < space.dof{1}.dvr_max % Right
        line ( [ hamilt.nip.max(1) hamilt.nip.max(1) ], ...
               [ space.dof{2}.dvr_min space.dof{2}.dvr_max ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.min(2) > space.dof{2}.dvr_min % Lower
        line ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ], ...
               [ hamilt.nip.min(2) hamilt.nip.min(2) ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.max(2) < space.dof{2}.dvr_max % Upper
        line ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max ], ...
               [ hamilt.nip.max(2) hamilt.nip.max(2) ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
end


%------------------------------------------------------------
% Densities projected on first coordinate (vertical plot)
%------------------------------------------------------------
function density_dvr1 ( step )       
global expect hamilt plots psi space
persistent rho_max

% Get projected densities where populations not too small
rho = cell (hamilt.coupling.n_eqs,1);
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        rho{m} = sum( abs(psi.dvr.grid_ND{m}).^2 .* space.dvr.weight_ND, 2 );
    end
end
        
% Get maximum density
% if step==1
    rho_max = 0;
    for  m=1:hamilt.coupling.n_eqs
        if expect.ind.pop{m}(step)>expect.min_pop
            rho_max = max ( rho_max, max(rho{m}) );
        end
    end
% end

% Plot projected densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        plot ( space.dvr.grid_1D{1},rho{m}/rho_max, ...
                 'LineStyle', '-', ...
                 'LineWidth', plots.style.line.thick, ...
                 'Color',     plots.style.colors(m,:)   );
    end
end

% Axes and labels
set ( gca, 'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
xlabel ( [ 'R_{', space.dof{1}.label, '}' ] )
ylabel ( [ '\rho (R_{', space.dof{1}.label, '})' ] )
if plots.density.range.on == false
    axis ( [ space.dof{1}.dvr_min space.dof{1}.dvr_max 0 1 ] )
else
    axis ( [ plots.density.range.x_min plots.density.range.x_max 0 1 ] )
end

% Absorbing boundary conditions
if ~isempty(hamilt.nip.grid_ND)
    if hamilt.nip.min(1) > space.dof{1}.dvr_min % Left
        line ( [ hamilt.nip.min(1) hamilt.nip.min(1) ], ...
               [ 0 1 ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.max(1) < space.dof{1}.dvr_max % Right
        line ( [ hamilt.nip.max(1) hamilt.nip.max(1) ], ...
               [ 0 1 ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end

end
%------------------------------------------------------------
% Densities projected on second coordinate (horizontal plot)
%------------------------------------------------------------
function density_dvr2 ( step )       
global expect hamilt plots psi space
persistent rho_max

rho  = cell (hamilt.coupling.n_eqs,1);

% Get projected densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        rho{m} = sum( abs(psi.dvr.grid_ND{m}).^2 .* space.dvr.weight_ND, 1 );
    end
end
        
% Get maximum density
% if step==1
rho_max = 0;
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        rho_max = max ( rho_max, max(rho{m}) );
    end
end
% end

% Plot projected densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        plot ( rho{m}/rho_max, space.dvr.grid_1D{2},...
            'LineStyle', '-', ...
            'LineWidth', plots.style.line.thick, ...
            'Color',     plots.style.colors(m,:)   );
    end
end

% Axes and labels
set ( gca, 'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
xlabel ( [ '\rho (R_{', space.dof{2}.label, '})' ] )
ylabel ( [ 'R_{', space.dof{2}.label, '}' ] )
if plots.density.range.on == false
    axis ( [ 0 1 space.dof{2}.dvr_min space.dof{2}.dvr_max ] )
else
    axis ( [ 0 1 plots.density.range.y_min plots.density.range.y_max ] )
end

% Absorbing boundary conditions
if ~isempty(hamilt.nip.grid_ND)
    if hamilt.nip.min(2) > space.dof{2}.dvr_min % Lower
        line ( [ 0 1 ], ...
               [ hamilt.nip.min(2) hamilt.nip.min(2) ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
    
    if hamilt.nip.max(2) < space.dof{2}.dvr_max % Upper
        line ( [ 0 1 ], ...
               [ hamilt.nip.max(2) hamilt.nip.max(2) ], ...
               'LineStyle', '--', ...
               'Color', 'k', ...
               'LineWidth', plots.style.line.thin)
    end
end

%------------------------------------------------------------
% Plot 2-dim momentum densities and kinetic energy surfaces
%------------------------------------------------------------
function density_fbr ( step )       
global expect hamilt info plots psi space

% Loop over individual wavefunctions
for  m=1:hamilt.coupling.n_eqs
    
    % Plot "mean trajectory" as solid curve in position space
    % Expectation value of R_1 vs. R_2
    if plots.density.expect.on
      plots.mask.ind{m} = find ( expect.ind.pop{m} > expect.min_pop );
      plot ( expect.ind.fbr{m}(plots.mask.ind{m},1), ...
           expect.ind.fbr{m}(plots.mask.ind{m},2), ...
          'LineStyle', '-', ...
          'LineWidth', plots.style.line.thick, ...
           'Color',     plots.style.colors(m,:) )
    end
       
    % Plot dotted contours of total enery functions
    % Creates total nonsense, maybe errors, for non-FFt grids, therefore restricted.
    % BAD: Uses internal object data.
    if plots.density.energy.on && isa(space.dof{1}, 'grid.fft') ...
       && isa(space.dof{2}, 'grid.fft')
        contour ( space.fbr.grid_ND{1} , ...
                  space.fbr.grid_ND{2}, ...
                  space.dof{1}.kin + space.dof{2}.kin, ...
                  linspace( plots.density.kin.min, plots.density.kin.max, plots.density.contour.nlev(2) ), ...
                  'LineStyle', ':', ...
                  'LineWidth', plots.style.line.thin, ...
                  'Color',     'k'  );
    end
    
    % If populations not too small, plot solid contours of densities
    if expect.ind.pop{m}(step)>expect.min_pop
        fbr = dvr2fbr(space.dof{1}, psi.dvr.grid_ND{m});
        fbr = dvr2fbr(space.dof{2}, fbr);
        contour ( space.fbr.grid_ND{1},    ...
                  space.fbr.grid_ND{2}, ...
                  abs(fbr).^2, ...
                  linspace(plots.density.rho_max.fbr/plots.density.contour.nlev(1), plots.density.rho_max.fbr, plots.density.contour.nlev(1)), ...   
                 'LineStyle', '-', ...
                 'LineWidth', plots.style.line.thin, ...
                 'Color',     plots.style.colors(m,:)   );
    end
    
end

axis ( [ space.fbr.grid_1D{1}(1) space.fbr.grid_1D{1}(end) ...
         space.fbr.grid_1D{2}(1) space.fbr.grid_1D{2}(end) ] )
set ( gca, 'XAxisLocation', 'top', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
title ( {info.header1;info.header2} )
if ~plots.density.marginals.on
    xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
end
ylabel ( [ 'P_{', space.dof{2}.label, '}' ] )

%------------------------------------------------------------
% Densities projected on first coordinate (vertical plot)
% Momentum space
%------------------------------------------------------------
function density_fbr1 ( step )       
global expect hamilt plots psi space
persistent rho_max

% Get projected densities where populations not too small
rho = cell (hamilt.coupling.n_eqs,1);
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        fbr = dvr2fbr(space.dof{1}, psi.dvr.grid_ND{m});
        fbr = dvr2fbr(space.dof{2}, fbr);
        rho{m} = sum( abs(fbr).^2, 2 );
    end
end
        
% Get maximum density
% if step==1
    rho_max = 0;
    for  m=1:hamilt.coupling.n_eqs
        if expect.ind.pop{m}(step)>expect.min_pop
            rho_max = max ( rho_max, max(rho{m}) );
        end
    end
% end

% Plot projected densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        plot ( space.fbr.grid_1D{1},rho{m}/rho_max, ...
                 'LineStyle', '-', ...
                 'LineWidth', plots.style.line.thick, ...
                 'Color',     plots.style.colors(m,:)   );
    end
end

% Axes and labels
set ( gca, 'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
xlabel ( [ 'P_{', space.dof{1}.label, '}' ] )
ylabel ( [ '\rho (P_{', space.dof{1}.label, '})' ] )
if plots.density.range.on == false
    axis ( [ space.dof{1}.fbr_min space.dof{1}.fbr_max 0 1 ] )
else
    axis ( [ plots.density.range.x_min plots.density.range.x_max 0 1 ] )
end

%------------------------------------------------------------
% Densities projected on second coordinate (horizontal plot)
% FBR plot
%------------------------------------------------------------
function density_fbr2 ( step )       
global expect hamilt plots psi space
persistent rho_max

rho  = cell (hamilt.coupling.n_eqs,1);

% Get projected densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        fbr = dvr2fbr(space.dof{1}, psi.dvr.grid_ND{m});
        fbr = dvr2fbr(space.dof{2}, fbr);
        rho{m} = sum( abs(fbr).^2, 1 );
    end
end
        
% Get maximum density
% if step==1
rho_max = 0;
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        rho_max = max ( rho_max, max(rho{m}) );
    end
end
% end

% Plot projected densities where populations not too small
for  m=1:hamilt.coupling.n_eqs
    if expect.ind.pop{m}(step)>expect.min_pop
        plot ( rho{m}/rho_max, space.fbr.grid_1D{2},...
            'LineStyle', '-', ...
            'LineWidth', plots.style.line.thick, ...
            'Color',     plots.style.colors(m,:)   );
    end
end

% Axes and labels
set ( gca, 'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'LineWidth',     plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',      plots.style.font.large, ...
           'FontWeight',    plots.style.font.heavy )
xlabel ( [ '\rho (P_{', space.dof{2}.label, '})' ] )
ylabel ( [ 'P_{', space.dof{2}.label, '}' ] )
if plots.density.range.on == false
    axis ( [ 0 1 space.dof{2}.fbr_min space.dof{2}.fbr_max ] )
else
    axis ( [ 0 1 plots.density.range.y_min plots.density.range.y_max ] )
end
