%----------------------------------------------------------
%
% Plot spectrum (Fourier transform of autocorrelation)
%
%----------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function spectrum 
global hamilt info plots time

% Plot total population (black curve)
if strcmp (info.program,'qm_bound')
    x = hamilt.eigen.eig_vals;
    y = ones(size(x));
    stem ( x,y, ...
        'LineStyle', '-', ...
        'LineWidth', plots.style.line.thick, ...
        'Color',     'k', ...
        'Marker', 'none')
elseif strcmp (info.program,'qm_propa')
    plot ( time.freq.grid, abs(time.spec.grid), ...
        'LineStyle', '-', ...
        'LineWidth', plots.style.line.thick, ...
        'Color',     'k' )
end
    
% Axes range
axis ( [ plots.density.tef.min plots.density.tef.max -0.1 1.1 ] )

% Labels
set ( gca, 'LineWidth',  plots.style.line.thick, ...
           'FontName',      plots.style.font.name,  ...
           'FontSize',   plots.style.font.large, ...
           'FontWeight', plots.style.font.heavy )   

xlabel ( 'energy' )
if info.program=='qm_propa'
    ylabel ( 'intensity' )
end

