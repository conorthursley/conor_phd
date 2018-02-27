%--------------------------------------------------------------------------
%
% Comparison of 'balanced truncation' and/or 'singular perturbation theory'
% and/or 'H2 optimal model reduction' versus dynamics in full dimensionality.
% For equations of motions as determined  by input variable "eom". #
% The input vectors "truncate", "singular" and/or "H2reduce" 
% contain a variable number of approximation degrees to be
% simulated; if one of those vectors is left empty on input, the 
% respective simulation method will be skipped. Graphical output will
% be displayed in figure specified by integer "fig".
% 
% switchoff=0: Do complete simulations
% switchoff>0: Skip calculation of bound states and matrix elements 
% switchoff>1: Skip calculation of A, B, N, C system matrices
% switchoff>2: Skip balancing transformation (BT method only)
%--------------------------------------------------------------------------


% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2013-16 Burkhard Schmidt, FU Berlin

function qm_BTversusH2 (eom, fig, truncate, singular, H2reduce, switchoff)

% Initializes general information and sets up log files.
init.info (mfilename('fullpath'));

% Initialize: Bound states, matrix elements, full dimensional, balancing
if switchoff<1 % Not necessary if energy/dipole matrices from external sources
    qm_setup();  qm_init; qm_bound  ;   qm_cleanup; 
    qm_setup();  qm_init; qm_matrix ;   qm_cleanup; 
end

if switchoff<2 % Not necessary if ABNCD from external sources
    qm_setup(); qm_init; qm_abncd   (eom); qm_cleanup;
end

% Control problem in full dimensionality as reference solution
qm_setup(); qm_init; qm_control (eom); qm_cleanup;
global bilinear control reduce info plots time
reference = control.y.forward;

ndim = size (control.x.forward,1);
for j =1:length (bilinear.label)
    mylegend {j,1} = strcat(int2str(ndim), ' (', bilinear.label{j}, ')');
end

% open figure
figure (fig); clf
info.program = 'qm_BTversusH2';
plot.logo;

%% Simple truncation
if ~isempty (truncate)
    subplot(1,3,1)
    hold on
    set ( gca, 'LineWidth',     plots.style.line.thick, ...
        'FontName',      plots.style.font.name,  ...
        'FontSize',      plots.style.font.large, ...
        'FontWeight',    plots.style.font.heavy )
    for len=1:control.y.dim
        plot (time.main.grid,reference(len,:)-reference(len,1),...
            'LineStyle', plots.style.patterns{1}, ...
            'LineWidth', plots.style.line.thin,...
            'Color',plots.style.colors(len,:))
    end
    
    % Balancing transformation
    if switchoff<3
        qm_setup; qm_init; qm_balance (eom);
        switchoff = 3;
    end
    
    %  Loop: different truncations
    for k=1:length(truncate)
        
        % Perform truncation and run in reduced dimensionality
        qm_truncate(eom,'t',truncate(k));
        qm_control (eom,'t',truncate(k));
        qm_cleanup;
        
        % Plot resulting populations
        figure(fig);
        subplot(1,3,1)
        global bilinear control reduce plots time
        for len=1:control.y.dim
            plot(time.main.grid,real(control.y.forward(len,:)-control.y.forward(len,1)),...
                'LineStyle', plots.style.patterns{k+1}, ...
                'LineWidth', plots.style.line.thin,...
                'Color',plots.style.colors(len,:))
        end
        
        for j =1:length (bilinear.label)
            mylegend {control.y.dim*k+j,1} = strcat(int2str(truncate(k)), ' (', bilinear.label{j}, ')');
        end
        
    end
    hold off
    legend (mylegend, 'Location','NorthEast', 'FontSize',plots.style.font.small)
    xlabel ('time t')
    ylabel ('observables y(t)-y(0)')
    switch lower(reduce.balance.A_stable)
        case 'ssu'
            title  ({'Balancing & Simple truncatation',['SSU stabilization: ' int2str(reduce.balance.A_split) ', B/N scaling: ' num2str(reduce.balance.BN_scale)],['ACF couple: ',int2str(reduce.balance.acf_couple),', Transform: ',upper(reduce.balance.transform)]})
        case 'evs'
            title  ({'Balancing & Simple truncatation',['EVS stabilization: ' num2str(reduce.balance.A_shift) ', B/N scaling: ' num2str(reduce.balance.BN_scale)],['Transform: ',upper(reduce.balance.transform)]})
    end
    drawnow
    
end

%% Singular perturbation
if ~isempty (singular)
    subplot(1,3,2)
    hold on
    set ( gca, 'LineWidth',     plots.style.line.thick, ...
        'FontName',      plots.style.font.name,  ...
        'FontSize',      plots.style.font.large, ...
        'FontWeight',    plots.style.font.heavy )
    for len=1:control.y.dim
        plot (time.main.grid,reference(len,:)-reference(len,1),...
            'LineStyle', plots.style.patterns{1}, ...
            'LineWidth', plots.style.line.thin,...
            'Color',plots.style.colors(len,:))
    end
    
    % Balancing transformation
    if switchoff<3
        qm_setup; qm_init; qm_balance (eom);
    end
    
    %  Loop: different truncations
    for k=1:length(singular)
        
        % Perform singular perturbation and run in reduced dimensionality
        qm_truncate(eom,'s',singular(k));
        qm_control (eom,'s',singular(k));
        qm_cleanup;
        
        % Plot resulting populations
        figure(fig);
        subplot(1,3,2)
        global bilinear control reduce plots time
        for len=1:control.y.dim
            plot(time.main.grid,real(control.y.forward(len,:)-control.y.forward(len,1)),...
                'LineStyle', plots.style.patterns{k+1}, ...
                'LineWidth', plots.style.line.thin,...
                'Color',plots.style.colors(len,:))
        end
        
        for j =1:length (bilinear.label)
            mylegend {control.y.dim*k+j,1} = strcat(int2str(singular(k)), ' (', bilinear.label{j}, ')');
        end
        
    end
    hold off
    legend (mylegend, 'Location','NorthEast', 'FontSize',plots.style.font.small)
    xlabel ('time t')
    ylabel ('observables y(t)-y(0)')
    switch lower(reduce.balance.A_stable)
        case 'ssu'
            title  ({'Balancing & Singular perturbation',['SSU stabilization: ' int2str(reduce.balance.A_split) ', B/N scaling: ' num2str(reduce.balance.BN_scale)],['ACF couple: ',int2str(reduce.balance.acf_couple),', Transform: ',upper(reduce.balance.transform)]})
        case 'evs'
            title  ({'Balancing & Singular perturbation',['EVS stabilization: ' num2str(reduce.balance.A_shift) ', B/N scaling: ' num2str(reduce.balance.BN_scale)],['Transform: ',upper(reduce.balance.transform)]})
    end
    drawnow
    
end

%% H2 optimal model reduction
if ~isempty (H2reduce)
    figure(fig)
    subplot(1,3,3)
    hold on
    set ( gca, 'LineWidth',     plots.style.line.thick, ...
        'FontName',      plots.style.font.name,  ...
        'FontSize',      plots.style.font.large, ...
        'FontWeight',    plots.style.font.heavy )
    
    for len=1:control.y.dim
        plot (time.main.grid,real(reference(len,:)-reference(len,1)),...
            'LineStyle', plots.style.patterns{1}, ...
            'LineWidth', plots.style.line.thin,...
            'Color',plots.style.colors(len,:))
    end
    
    for k=1:length(H2reduce)
        
        % Perform dimension reduction and run in reduced dimensionality
        qm_setup();
        qm_init;
        qm_H2model (eom,    H2reduce(k));
        qm_control (eom,'h',H2reduce(k));
        qm_cleanup;
        
        % Plot resulting populations
        figure(fig);
        subplot(1,3,3)
        global bilinear control reduce plots time
        for len=1:control.y.dim
            plot(time.main.grid, real(control.y.forward(len,:)-control.y.forward(len,1)),...
                'LineStyle', plots.style.patterns{k+1}, ...
                'LineWidth', plots.style.line.thin,...
                'Color',plots.style.colors(len,:))
        end
        
        for j =1:length (bilinear.label)
            mylegend {control.y.dim*k+j,1} = strcat(int2str(H2reduce(k)), ' (', bilinear.label{j}, ')');
        end
        
    end
    hold off
    legend (mylegend, 'Location','NorthEast', 'FontSize',plots.style.font.small)
    xlabel ('time t')
    ylabel ('observables y(t)-y(0)')
    switch lower(reduce.H2model.A_stable)
        case 'ssu'
            title  ({'H2 model reduction',['SSU stabilization: ' int2str(reduce.H2model.A_split) ', B/N scaling: ' num2str(reduce.H2model.BN_scale)],['Tolerance: ',num2str(reduce.H2model.conv_tol),', Max. iterations: ',int2str(reduce.H2model.max_iter)]})
        case 'evs'
            title  ({'H2 model reduction',['EVS stabilization: ' num2str(reduce.H2model.A_shift) ', B/N scaling: ' num2str(reduce.H2model.BN_scale)],['Tolerance: ',num2str(reduce.H2model.conv_tol),', Max. iterations: ',int2str(reduce.H2model.max_iter)]})
    end
    
end

% Output clock/date/time
util.clock;

% Save figure
saveas (gca, int2str(fig), 'fig')

end


