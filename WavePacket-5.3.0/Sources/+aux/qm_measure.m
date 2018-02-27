%--------------------------------------------------------------------------
%
% Reads the solution x(t) from bilinear control equation and
% undoes balanced truncation transformation and
% reconstructs the density operator (in matrix form) and
% calculates various correlation measures which may be of
% interest for quantum information theory
%
% purity   : tr(rho^2)
% impurity : 1 - tr(rho^2)
% entropy  : -tr(rho*log(rho))    (von Neumann)
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2013 Burkhard Schmidt, Jeremy Rodriguez, Ulf Lorenz
%
% see the README file for license details.

% Input variable "filename" typically specifying physical context, e.g.
% 'lvne' for Liouville-von-Neumann
% 'tdse' for time-dependent Schroedinger equation
% Input variable "filenumber" specifying the balancing/truncation scheme
% =0 - unbalanced, untruncated
% =1 - balanced untruncated
% >1 - truncated to this many equations 

function qm_measure(filename, filenumber)

% Main variables are global throughout;
global control correlate plots time 

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

util.disp (' ')
util.disp ('-------------------------------------------------------------')
util.disp (' Reads the solution x_t from bilinear control equation and   ')
util.disp (' undoes balanced truncation transformation and               ')
util.disp (' reconstructs the density operator (in matrix form) and      ')
util.disp (' calculates various correlation measures                     ')
util.disp ('                                                             ')
util.disp ('-------------------------------------------------------------')
util.disp (' ')

%% Initialize everything

% Initialize temporal discretization
init.timesteps;

% Load A, B, N, and C, D matrices, initial/equilib state/density, etc.
load ([filename '_' int2str(filenumber)]);

% Load time dependence of x-vector
load (['uxy_' filename '_' int2str(filenumber)])

% Prepare figure legend
correlate.legend = {};

% Prepare calculation of purity/impurity
if isfield(control.correlation, 'purity')
    correlate.purity = zeros(time.main.n, 1);
    correlate.impure = zeros(time.main.n, 1);
    correlate.legend = [correlate.legend 'purity = tr(\rho^2)'];
    correlate.legend = [correlate.legend 'impurity = 1-purity'];
end

% Prepare calculation of von Neumann entropy
if isfield(control.correlation, 'entropy')
    correlate.entropy = zeros(time.main.n, 1);
    correlate.max_ent = zeros(time.main.n, 1);
    correlate.legend = [correlate.legend 'entropy = -tr(\rho log \rho)'];
    correlate.legend = [correlate.legend 'max(entropy) = log(N)'];
end

% Prepare calculation of selected coherences
if isfield(control.correlation, 'coherence')
            correlate.coherence = cell(length(control.correlation.coherence), 1);
            correlate.max_coher = cell(length(control.correlation.coherence), 1);
            
            for len = 1:length(control.correlation.coherence)
                ii = control.correlation.coherence{len}(1)+1;
                jj = control.correlation.coherence{len}(2)+1;
                
                correlate.coherence{len} = zeros(time.main.n, 1);
                correlate.max_coher{len} = zeros(time.main.n, 1);
                
                correlate.legend=[correlate.legend ['coherence |' num2str(ii-1) '><' num2str(jj-1) '|']];
                correlate.legend=[correlate.legend  'maximum thereof'];
            end
    
end

%% Main loop over time steps
for step = 1:time.main.n
    
    % Re-construct density, i.e. undo balancing/truncation
    rho = oct.reconstruct (x.t(:,step)+x.e, filename, filenumber);
    
    % Calculate purity / impurity
    if  isfield(control.correlation, 'purity')
        correlate.purity(step) = real(trace(rho^2));
        correlate.impure(step) = 1 - correlate.purity(step);
    end
    
    % Calculate von Neumann entropy
    if isfield(control.correlation, 'entropy')        
        lambda = real(eig(rho));
        S = 0;
        for ii=1:length(lambda)
           if lambda(ii)>0 % avoid zero eigenvalues
               S = S - lambda(ii)*log(lambda(ii));
           end
        end
        correlate.entropy(step) = S; 
        correlate.max_ent(step) = log(size(rho,2));
    end
    
    % Calculate selected coherences
    if isfield(control.correlation, 'coherence')
        for len = 1:length(control.correlation.coherence)
            ii = control.correlation.coherence{len}(1)+1;
            jj = control.correlation.coherence{len}(2)+1;
            
            % Take absolute value to suppress fast oscillations
            correlate.coherence{len}(step) = abs(rho(ii,jj));
            
            % Get upper bound from Schwartz inequality
            correlate.max_coher{len}(step) = sqrt(abs(rho(ii,ii)*rho(jj,jj)));
        end
    end
    
end

%% Graphical representations of the results
figure(21);
clf;
plot.logo
hold on

% Plot purity/impurity versus time
if isfield(control.correlation, 'purity')
    plot (time.main.grid, correlate.purity,...
        'LineStyle', '-', ...
        'LineWidth', plots.style.line.thick,...
        'Color',     plots.style.colors(1,:))
    plot (time.main.grid, correlate.impure,...
        'LineStyle', ':', ...
        'LineWidth', plots.style.line.thick,...
        'Color',     plots.style.colors(1,:))
end

% Plot von Neumann entropy versus time
if isfield(control.correlation, 'entropy')
    plot (time.main.grid, correlate.entropy,...
        'LineStyle', '-', ...
        'LineWidth', plots.style.line.thick,...
        'Color',     plots.style.colors(2,:))
    plot (time.main.grid, correlate.max_ent,...
        'LineStyle', ':', ...
        'LineWidth', plots.style.line.thick,...
        'Color',     plots.style.colors(2,:))
end

% Plot selected coherences versus time
if isfield(control.correlation, 'coherence')
    for len = 1:length(control.correlation.coherence)
        plot (time.main.grid, correlate.coherence{len},...
            'LineStyle', '-', ...
            'LineWidth', plots.style.line.thick,...
            'Color',     plots.style.colors(2+len,:))
        plot (time.main.grid, correlate.max_coher{len},...
            'LineStyle', ':', ...
            'LineWidth', plots.style.line.thick,...
            'Color',     plots.style.colors(2+len,:))
    end
end

hold off
set ( gca, 'LineWidth',     plots.style.line.thick, ...
    'FontName',      plots.style.font.name,  ...
    'FontSize',      plots.style.font.large, ...
    'FontWeight',    plots.style.font.heavy )

legend (correlate.legend)
xlabel('time')
ylabel('correlation measures')
title([ matrix.title int2str(size(x.t,1)) ' coupled ODEs'])


