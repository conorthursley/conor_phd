%------------------------------------------------------------------------------
%
% Initial settings of plot parameters
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2011 Ulf Lorenz
%               2009 Burkhard Schmidt
%               2012 Ulf Lorenz
%
% see the README file for license details.

function plots
global plots


%%  General appearance of plots 

screen = get (0, 'ScreenSize');          % Screen size
ScreenHeight = screen(4);                % Don't use width for double screens

% If matlab is started from the command line (matlab -nodisplay), the
% above code screws up. In this case, supply a standard value for the
% dimensions. 
if ScreenHeight < 100
	ScreenHeight = 768;
end

% RGB colors for individual SEs: Matlab's color scheme introduced in R2014a
plots.style.colors = [...
         0    0.4470    0.7410; ...
    0.8500    0.3250    0.0980; ...
    0.9290    0.6940    0.1250; ...
    0.4940    0.1840    0.5560; ...
    0.4660    0.6740    0.1880; ...
    0.3010    0.7450    0.9330; ...
    0.6350    0.0780    0.1840];
plots.style.colors = repmat(plots.style.colors,10,1); % repeat 10 times

% Patterns for curves/lines
plots.style.patterns = {'-'; ...         % 1: solid (default) 
                        '--';  ...       % 2: dashed
                        ':'; ...         % 3: dotted
                        '-.'; ...        % 4: dash-dotted
                        '-'; ...         % 5: solid (default)
                        '--';  ...       % 6: dashed
                        ':'; ...         % 7: dotted
                        '-.'; ...        % 8: dash-dotted
                        '-'; ...         % 9: solid (default)
                        '--';  ...       % 10: dashed
                        ':'; ...         % 11: dotted
                        '-.'};           % 12: dash-dotted

 
% Line widths
plots.style.line.thin  = 1;              % energy contours/curves
plots.style.line.thick = 2;              % axes, density contours/curves
plots.style.line.extrathick = 4;         % polar plots

% Font sizes and weights
plots.style.font.large = round(ScreenHeight/125);
plots.style.font.small = round(ScreenHeight/125)-2;
plots.style.font.heavy = 'bold';         
plots.style.font.light = 'light';        
plots.style.font.name  = 'Helvetica';

% Marker sizes 
plots.style.marker.large = round(ScreenHeight/400);
plots.style.marker.small = round(ScreenHeight/400)-2;

%% Plot densities evolving in time

% Factor by which densities are divided before displaying them
plots.density.rho_max.dvr = [];
plots.density.rho_max.fbr = [];

% Minimum and maximum ranges to plot potentials
plots.density.pot.min = [];
plots.density.pot.max = [];

% Toggle plot
plots.density.on = true;   

% Position and size of plot window: Format 9:9
plots.density.size.left   = round(  ScreenHeight/16);     % Left border
plots.density.size.lower  = round(  ScreenHeight/16);     % Lower border
plots.density.size.width  = round(9*ScreenHeight/16);     % Width
plots.density.size.height = round(9*ScreenHeight/16);     % Height

plots.density.size.width =  4*round(plots.density.size.width/4);
plots.density.size.height = 4*round(plots.density.size.height/4);

% Appearence of density plots
plots.density.type          = 'contour'; % 'contour' draws a contour plot of WF/Wigner trafo
                                         % 'surface' does the same with a rendered surface
                                         % 'reduced' produces reduced Wigner transforms
                                         % 'curve' just draws the 1D wavefunction
                                         % 'polar' makes a polar plot of the 1D wavefunction
plots.density.hold          = false;     % Hold previous densities or overwrite them
plots.density.representation= 'dvr';     % Position (dvr) or momentum (fbr) densities
plots.density.curve.what    = 'density'; % Choice of what to show: realwav|imagwav|density
plots.density.contour.nlev  = [30 15];   % Number of contours: density/energy
plots.density.contour.min = [];          % Lower contour line for density plots
plots.density.contour.max = [];          % Upper contour line for density plots
plots.density.surface.view  = [60 75];   % View point for surface plot: [az el]
plots.density.surface.look  = [ 1  1];   % Look of surface plot [shading lighting]
plots.density.surface.lite  = [45 45];   % Angle of light for surface plot: [az el]
plots.density.surface.color.on  = false; % manual setting of color ranges
plots.density.surface.color.min = [];    % minimum energy value passed to caxis()
plots.density.surface.color.max = [];    % maximum energy value passed to caxis()
plots.density.surface.color.map = 'default';    % color map
plots.density.expect.on     = false;     % Toggle expectation values 
plots.density.energy.on     = true;      % Toggle energy surfaces
plots.density.marginals.on  = true;      % Toggle marginal function
plots.density.logo.on       = true;      % Logos in all four corners
plots.density.kin.min = [];              % Minimum value when plotting kin. energy
plots.density.kin.max = [];              % Maximum value when plotting kin. energy
plots.density.pot.min = [];              % Minimum value when plotting pot. energy
plots.density.pot.max = [];              % Maximum value when plotting pot. energy
plots.density.range.on      = false;     % Toggle manual setting of plot range
plots.density.range.x_min   = [];        % Minimum of first coordinate
plots.density.range.x_max   = [];        % Maximum of first coordinate
plots.density.range.y_min   = [];        % Minimum of second/momentum coordinate
plots.density.range.y_max   = [];        % Maximum of second/momentum coordinate

% Export to video file: Currently only mp4-format supported
plots.density.export.on     = true;      % Toggle video export
plots.density.export.file   = [];        % Set the output to a custom filename (w/o suffix)
plots.density.export.images = false;     % Export movie as a series of images


%% Plot expectation values evolving in time
      
% Toggle plot
plots.expect.on = true;   

% Position and size of plot: Format 7:9
plots.expect.size.left   = round(11*ScreenHeight/16);     % Left border
plots.expect.size.lower  = round(01*ScreenHeight/16);     % Lower border
plots.expect.size.width  = round(07*ScreenHeight/16);     % Width
plots.expect.size.height = round(09*ScreenHeight/16);     % Height

% Appearence of plot
plots.expect.errorbar.on    = false;     % Toggle plotting errorbars
plots.expect.population.min = -0.1;      % Minimum of norm plot
plots.expect.population.max = 1.1;       % Maximum of norm plot
plots.expect.energies.min   = [];        % Lower bound of the energy plot
plots.expect.energies.max   = [];        % Upper bound of the energy plot
plots.expect.legends.on     = true;      % Toggle legends
plots.expect.logo.on        = true;      % Logos in all four corners

% Export to graphics file: jpg/tiff/epsc/pdf, etc ...
plots.expect.export.on      = false;     % Toggle graphics export    
plots.expect.export.file    = [];        % Set custom filename (suffix determines image type)


%% Plot spectrum (Fourier transform of autocorrelation)
      
% Toggle plot
plots.spectrum.on = false;   

% Position and size of plot 
plots.spectrum.size.left   = round(01  *ScreenHeight/16);     % Left border
plots.spectrum.size.lower  = round(11.5*ScreenHeight/16);     % Lower border
plots.spectrum.size.width  = round(16  *ScreenHeight/16);     % Width
plots.spectrum.size.height = round(03  *ScreenHeight/16);     % Height

% Appearence of plot
plots.spectrum.logo.on     = true;       % Logos in all four corners

% Export to graphics file: jpg/tiff/epsc/pdf, etc ...
plots.spectrum.export.on   = false;      % Toggle graphics export    
plots.spectrum.export.file = [];         % Set custom filename (suffix determines image type)


%% Plot of control input/output/observables
plots.control.size.left   = round(   ScreenHeight/16);     % Left border
plots.control.size.lower  = round(   ScreenHeight/16);     % Lower border
plots.control.size.width  = round( 7*ScreenHeight/16);     % Width
plots.control.size.height = round(12*ScreenHeight/16);     % Height
