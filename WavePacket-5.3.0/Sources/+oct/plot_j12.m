% Plot three control functionals during iteration
% (1) Target
% (2) Costs
% (3) Total
%
%
% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015-16 Burkhard Schmidt
%
% see the README file for license details.

function plot_j12 ( iter )

global plots control

% Initialize figure; draw logos, and resize
if iter==1
    h = figure(8);
    clf;
    plot.logo;
    set(h,'units','pixels', ...
        'position',[ ...
        plots.control.size.left+plots.control.size.width*1.1 ...
        plots.control.size.lower ...
        plots.control.size.width ...
        plots.control.size.height] );
else
    figure(8)
end

% Top: target functional
subplot(3,1,1)
plot(1:iter,control.j.target(1:iter),'o', ...
    'LineStyle', '-', ...
    'LineWidth', plots.style.line.thin)
set ( gca, 'LineWidth',     plots.style.line.thick, ...
    'FontName',      plots.style.font.name,  ...
    'FontSize',      plots.style.font.large, ...
    'FontWeight',    plots.style.font.heavy )
xlabel('iterations')
ylabel('target J_1')
title (control.title);

% Middle: Cost functional
subplot(3,1,2)
plot(1:iter,control.j.cost(1:iter),'o', ...
    'LineStyle', '-', ...
    'LineWidth', plots.style.line.thin)
set ( gca, 'LineWidth',     plots.style.line.thick, ...
    'FontName',      plots.style.font.name,  ...
    'FontSize',      plots.style.font.large, ...
    'FontWeight',    plots.style.font.heavy )
xlabel('iterations')
ylabel('cost J_2')

% Bottom: Total functional
subplot(3,1,3)
plot(1:iter,control.j.total(1:iter),'o', ...
    'LineStyle', '-', ...
    'LineWidth', plots.style.line.thin)
hold on
plot(1:iter,control.j.addup(1:iter),'o', ...
    'LineStyle', ':', ...
    'LineWidth', plots.style.line.thin)
hold off
xlabel('iterations')
ylabel('optimize total')
set ( gca, 'LineWidth',     plots.style.line.thick, ...
    'FontName',      plots.style.font.name,  ...
    'FontSize',      plots.style.font.large, ...
    'FontWeight',    plots.style.font.heavy )

% Update figure window
drawnow

saveas (8,'j12.fig')
saveas (8,'j12.jpg')

