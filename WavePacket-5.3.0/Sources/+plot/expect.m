%--------------------------------------------------------------------------
%
% Compose animated figure from several subplots
% displaying the time evolution of various expectation values
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2009-2016 Burkhard Schmidt
%
% see the README file for license details.

function expect (wide)
global expect hamilt info plots space time
   
if wide % Wide format: 16:9
    w=16; h=09; o=9/16;
else % Narrow format: 7:9
    w=07; h=09; o=0;
end

% Find steps where individual populations are non-negligible
for m=1:hamilt.coupling.n_eqs
    plots.mask.ind{m} = find ( expect.ind.pop{m} > expect.min_pop );
end

% Fake: Find steps already existing
plots.mask.tot = find ( expect.tot.pop > expect.min_pop );

% TDSE = Time dependent Schrödinger equation
if strcmpi (info.program,'qm_propa')
    
    % Populations/projections vs. time
    subplot ( 'Position', [o+1/w (1+14/3)/h 5/w 7/(3*h)] )
    cla;
    population
    if wide; set ( gca, 'YAxisLocation', 'right'); end
    
    % Expectation values: Energies versus time
    subplot ( 'Position', [o+1/w (1+7/3)/h 5/w 7/(3*h)] )
    cla;
    energies
    if wide; set ( gca, 'YAxisLocation', 'right'); end
    
    if time.efield.n_pulse == 0
        
        % Autocorrelation function versus time
        subplot ( 'Position', [o+1/w 1/h 5/w 7/(3*h)] )
        cla;
        correlation
        if wide; set ( gca, 'YAxisLocation', 'right'); end
        
    else
        
        % External electric field versus time
        subplot ( 'Position',[o+1/w 1/h 5/w 7/(3*h)] )
        cla;
        efield
        if wide; set ( gca, 'YAxisLocation', 'right'); end
        
    end
    

% TISE = Time independent Schrödinger equation
elseif strcmpi (info.program,'qm_bound')
    
    % Populations/projections and energies
    if isfield(space, 'amo')
        subplot ( 'Position', [o+1/w 5/h 5/w 3/h] )
        population
        if wide; set ( gca, 'YAxisLocation', 'right'); end
        
        subplot ( 'Position', [o+1/w 1/h 5/w 3/h] )
        energies
        if wide; set ( gca, 'YAxisLocation', 'right'); end
    else
    
    % Energies only
        subplot ( 'Position', [o+1/w 1/h 5/w 7/h] )
        energies
        if wide; set ( gca, 'YAxisLocation', 'right'); end
    end

end
    
%--------------------------------------------------------------------
% Plot populations/projections vs. time
%--------------------------------------------------------------------
function population
global expect hamilt info plots space time

% Plot total population (black curves, solid)
plot ( time.main.grid ( plots.mask.tot ), ...
    expect.tot.pop ( plots.mask.tot ), ...
    'LineStyle', plots.style.patterns{1}, ...
    'LineWidth', plots.style.line.thick, ...
    'Color',     'k' )
hold on

% Plot additional multiplicative operators (black curves, different linestyles)
if isfield(space, 'amo')
    for p = 1:length(space.amo)
        plot ( time.main.grid ( plots.mask.tot ), ...
            expect.tot.amo ( plots.mask.tot, p ), ...
            'LineStyle', plots.style.patterns{p+1}, ...
            'LineWidth', plots.style.line.thick, ...
            'Color',     'k' )
    end
end

if plots.expect.legends.on==1
    labels{1} = 'Population';
    if isfield(space, 'amo')
        for p = 1:length(space.amo)
            labels{end+1} = space.amo{p}.label;
        end
    end
end

% Plot populations for individual wavefunctions (colored curves)
n = hamilt.coupling.n_eqs;
if n>1
    
    % Plot populations
    for m=1:n
        plot ( time.main.grid    ( plots.mask.ind{m}), ...
               expect.ind.pop{m} ( plots.mask.ind{m}), ...
               'LineStyle', plots.style.patterns{1}, ...
               'LineWidth', plots.style.line.thick, ...
               'Color',     plots.style.colors(m,:))
    end
    
    % Legend explaining the colors
    if plots.expect.legends.on==1
        for m=1:n
            if strcmpi (hamilt.coupling.representation,'dia')
                labels{end+1} = hamilt.coupling.labels{m};
            elseif strcmpi (hamilt.coupling.representation,'adi')
                labels{end+1} = int2str(m);
            end
        end
    end
    
    % Plot projections
    if isfield(space, 'amo')
        for p = 1:length(space.amo)  
            for m=1:n
                plot ( time.main.grid    ( plots.mask.ind{m}), ...
                    expect.ind.amo{m} ( plots.mask.ind{m}, p), ...
                    'LineStyle', plots.style.patterns{p+1}, ...
                    'LineWidth', plots.style.line.thick, ...
                    'Color',     plots.style.colors(m,:))
            end
        end
    end

end

% Draw legend
legend(labels,'Location','SouthWest')

% Axes and labels
axis ( [ 0 time.main.total plots.expect.population.min plots.expect.population.max ] );
set ( gca, 'XAxisLocation', 'top', ...
           'LineWidth',  plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',   plots.style.font.large, ...
           'FontWeight', plots.style.font.heavy )   

if strcmpi (info.program,'qm_propa') % TDSE: Loop over time steps
     xlabel ('t')
elseif strcmpi (info.program,'qm_bound') % TISE: Loop over eigenstates
     xlabel ('n')
end

if hamilt.coupling.n_eqs==1 
    ylabel ( 'N(t)' )
else
    if strcmpi (hamilt.coupling.representation,'adi')
        ylabel ( 'N_{adi}(t)' )
    elseif strcmpi (hamilt.coupling.representation,'dia')
        ylabel ( 'N_{dia}(t)' )
    end
end

hold off

%------------------------------------------------------------------------------------
% Plot expectation values: sum of energies versus time
%------------------------------------------------------------------------------------
function energies
global expect hamilt info plots time uncert

%% Plot energies for total wavefunction (black curves)

% Total energy
plot ( time.main.grid(plots.mask.tot), ...
       expect.tot.tot(plots.mask.tot), ...
       'LineStyle','-', ...
       'LineWidth',plots.style.line.thick, ... 
       'Color','k')
hold on

% Potential energy
if hamilt.coupling.n_eqs==1
    if plots.expect.errorbar.on
        errorbar ( time.main.grid(plots.mask.tot), ...
            expect.tot.pot   (plots.mask.tot), ...
            uncert.ind.pot{1}(plots.mask.tot), ...  % ind{1}=tot for n_eqs=1
            'LineStyle','--', ...
            'LineWidth',plots.style.line.thick, ...
            'Color','k')
    else
        plot ( time.main.grid(plots.mask.tot), ...
            expect.tot.pot   (plots.mask.tot), ...
            'LineStyle','--', ...
            'LineWidth',plots.style.line.thick, ...
            'Color','k')
    end
else
    plot ( time.main.grid(plots.mask.tot), ...
        expect.tot.pot(plots.mask.tot), ...
        'LineStyle','--', ...
        'LineWidth',plots.style.line.thick, ...
        'Color','k')
end

% Kinetic energy
if hamilt.coupling.n_eqs==1
    if plots.expect.errorbar.on
        errorbar ( time.main.grid(plots.mask.tot), ...
            expect.tot.kin   (plots.mask.tot), ...
            uncert.ind.kin{1}(plots.mask.tot), ... % ind{1}=tot for n_eqs=1
            'LineStyle',':', ...
            'LineWidth',plots.style.line.thick, ...
            'Color','k')
    else
        plot ( time.main.grid(plots.mask.tot), ...
            expect.tot.kin   (plots.mask.tot), ...
            'LineStyle',':', ...
            'LineWidth',plots.style.line.thick, ...
            'Color','k')
    end
else
    plot ( time.main.grid(plots.mask.tot), ...
        expect.tot.kin(plots.mask.tot), ...
        'LineStyle',':', ...
        'LineWidth',plots.style.line.thick, ...
        'Color','k')
end
   
% Legend explaining the line styles
if plots.expect.legends.on==true
    legend('<E>','<V>','<T>','Location','SouthWest')
end
    
%% Plot all/potential/kinetic energies for individual wavefunctions (colored curves)
if hamilt.coupling.n_eqs>1
    for m=1:hamilt.coupling.n_eqs
        plot ( time.main.grid   (plots.mask.ind{m}), ...
            expect.ind.all{m}(plots.mask.ind{m}), ...
            'LineStyle','-', ...
            'LineWidth',plots.style.line.thick, ...
            'Color',plots.style.colors(m,:))
        if plots.expect.errorbar.on
            errorbar ( time.main.grid   (plots.mask.ind{m}), ...
                expect.ind.pot{m}(plots.mask.ind{m}), ...
                uncert.ind.pot{m}(plots.mask.ind{m}), ...
                'LineStyle','--', ...
                'LineWidth',plots.style.line.thick, ...
                'Color',plots.style.colors(m,:))
            errorbar ( time.main.grid   (plots.mask.ind{m}), ...
                expect.ind.kin{m}(plots.mask.ind{m}), ...
                uncert.ind.kin{m}(plots.mask.ind{m}), ...
                'LineStyle',':', ...
                'LineWidth',plots.style.line.thick, ...
                'Color',plots.style.colors(m,:))
        else
            plot ( time.main.grid   (plots.mask.ind{m}), ...
                expect.ind.pot{m}(plots.mask.ind{m}), ...
                'LineStyle','--', ...
                'LineWidth',plots.style.line.thick, ...
                'Color',plots.style.colors(m,:))
            plot ( time.main.grid   (plots.mask.ind{m}), ...
                expect.ind.kin{m}(plots.mask.ind{m}), ...
                'LineStyle',':', ...
                'LineWidth',plots.style.line.thick, ...
                'Color',plots.style.colors(m,:))
        end
    end
end

%% Axes, labels, etc
lower = plots.density.pot.min;
upper = max(plots.density.tef.max, plots.density.kin.max);
if ~isempty(plots.expect.energies.min)
    lower = plots.expect.energies.min;
end
if ~isempty(plots.expect.energies.max)
    upper = plots.expect.energies.max;
end
axis ( [ 0 time.main.total lower upper ] )
set ( gca, 'LineWidth',  plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',   plots.style.font.large, ...
           'FontWeight', plots.style.font.heavy )

if strcmpi (info.program,'qm_propa') % TDSE: Loop over time steps
     set ( gca, 'XTick',      [-123456789 123456789]) % Fake: suppress tick labels
elseif strcmpi (info.program,'qm_bound') % TISE: Loop over eigenstates
     xlabel ('n')
end
if hamilt.coupling.n_eqs==1 
    ylabel ( '<V,T,E>' )
else
    if strcmpi (hamilt.coupling.representation,'adi')
        ylabel ( '<V>, <T>, <E>_{adi}' )
    elseif strcmpi (hamilt.coupling.representation,'dia')
        ylabel ( '<V>, <T>, <E>_{dia}' )
    end
end
hold off

%--------------------------------------------------------------------
% Plot autocorrelation vs. time
%--------------------------------------------------------------------
function correlation
global plots time

% Not too many time steps: plot colored curve
if time.sub.n * time.main.n<=1000

    % Get modulus of autocorrelation function
    rho = abs   ( time.acf.grid ) .^2;

    % Get phase of autocorrelation function: Map interval [-pi,pi] into [0,1]
    phi = angle ( time.acf.grid ) / (2*(pi+0.001)) + 1/2;

    % Plot autocorrelation function (colored curve)
    plot.color ( time.sub.grid, rho, phi, [1 1 0]/2, plots.style.line.thick, 0, 0 )

% Too many time steps: plot black curve
else
    plot ( time.sub.grid, abs(time.acf.grid), ...
        'LineStyle', '-', ...
        'LineWidth', plots.style.line.thick, ...
        'Color',     'k' )

end

% Axes and labels
axis ( [ 0 time.main.total -0.1 1.1 ] )
set ( gca, 'LineWidth',  plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',   plots.style.font.large, ...
           'FontWeight', plots.style.font.heavy )   
xlabel ('t')
ylabel ('C(t)' )

%--------------------------------------------------------------------
% Plot external electric field vs. time
%--------------------------------------------------------------------
function efield
global plots time

if time.efield.dressed
    factor = 2;
else
    factor = 1;
end

% Plot x-component, y-component, magnitude of electric field
if any(time.efield.polar~=pi/2)
    plot ( time.sub.grid, factor * real(time.efield.grid.x), ...
        'LineStyle', '--', ...
        'LineWidth', plots.style.line.thick, ...
        'Color',     'k' )
    hold on
end
if any(time.efield.polar~=0)
    plot ( time.sub.grid, factor * real(time.efield.grid.y), ...
        'LineStyle', ':', ...
        'LineWidth', plots.style.line.thick, ...
        'Color',     'k' )
    hold off
end

% Legend explaining the line styles
if plots.expect.legends.on==1 && any(time.efield.polar~=pi/2) && any(time.efield.polar~=0)
    legend('F_x','F_y','Location','SouthWest')
end

% Axes and labels
if time.efield.dressed
    axis ( [ 0 time.main.total          -0.1             +max(time.efield.ampli)*1.1 ] )
else
    line ([0 time.main.total],[0 0],'Color','k','LineWidth', plots.style.line.thin)
    axis ( [ 0 time.main.total -max(time.efield.ampli)*1.1 +max(time.efield.ampli)*1.1 ] )
end
set ( gca, 'LineWidth',  plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',   plots.style.font.large, ...
           'FontWeight', plots.style.font.heavy )   
xlabel ('t')
ylabel ('F(t)' )
