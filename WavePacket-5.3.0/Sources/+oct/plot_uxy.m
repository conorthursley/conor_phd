% Plot evolution of input/state/output of bilinear control problem 
% in one figure with 2 or 3 subplots 
% (1) input=control field(s) u(t) - not drawn for the field-free case - 
% (2) state vector x(t)
% (3) output=observable(s) y(t) 
%
% Depending on variable 'action' the following is achieved:
% 'open'     : open figure/subplots and prepare legends
% 'equilib'  : draw equilibrium values (dashed lines)
% 'forward'  : draw forward propagation (solid lines)
% 'backward' : draw backward propagation (dotted lines)
% 'clear'    : clear figure/subplots 
%
% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2015-16 Burkhard Schmidt
%
% see the README file for license details.

function plot_uxy ( action, first, last )

global bilinear plots control time
persistent numPlots writerObj
control.x.max = 7;

switch lower(action)
    
    case 'open'
        
        % Get the number of subplots: 2 (without control) or 3 (with control)
        numPlots = 2 + isfield(control,'u');
        
        % Initialize figure and resize
        h = figure(7);
        clf;
        set(h,'units','pixels', ...
            'position',[ ...
            plots.control.size.left ...
            plots.control.size.lower ...
            plots.control.size.width ...
            plots.control.size.height] );
        
        % Optionally open a movie file
        if control.plot.mov
            util.disp ('Opening animated uxy plot file : uxy.mp4')
            writerObj = VideoWriter ('uxy', 'MPEG-4');
            open(writerObj);
        end
        
        % Prepare legends for u, x, y
        if isfield(control,'u')
            control.u.legends = cell(control.u.dim,1);
            for i=1:control.u.dim
                control.u.legends{i}=['u_{' int2str(i) '}(t)'];
            end
        end
        
        control.x.legends = cell(min(control.x.dim,control.x.max),1);
        for i=1:min(control.x.dim,control.x.max)
            control.x.legends{i}=['x_{' int2str(i) '}(t)'];
        end
        
        control.y.legends = bilinear.label;

        %% Initialize figure with 2 or 3 subplots
    case ('equilib')
        figure(7);
        
        % Clear figure and draw logos in all four corners
        clf;
        plot.logo;

        % (1) Input (=control fields) u(t)
        if isfield(control,'u')
            subplot(numPlots,1,1);
            
            axis ( [[control.t.steps(1) control.t.steps(end)] [-1 1]*2*abs(time.efield.ampli) ] )
            hold on
            for len=1:control.u.dim
                line([control.t.steps(1) control.t.steps(end)],[0 0],...
                    'LineStyle', ':', ...
                    'LineWidth', plots.style.line.thin,...
                    'Color',plots.style.colors(len,:))
            end
            set ( gca, 'LineWidth',     plots.style.line.thick, ...
                'FontName',      plots.style.font.name,  ...
                'FontSize',      plots.style.font.large, ...
                'FontWeight',    plots.style.font.heavy )
            xlabel('time')
            ylabel('control field(s) u(t)')
            box on
            legend(control.u.legends,'Location','NorthWest')
            
        end
        
        % (2) State vector x(t)
        subplot(numPlots,1,numPlots-1);
        
        axis ( [ [control.t.steps(1) control.t.steps(end)] [-0.1 1.1]] )
        hold on
        for len=1:min(control.x.dim,control.x.max)
            line([control.t.steps(1) control.t.steps(end)],real(control.x.equilib(len))*[1 1],...
                'LineStyle', '--', ...
                'LineWidth', plots.style.line.thin,...
                'Color',plots.style.colors(len,:))
        end
        set ( gca, 'LineWidth',     plots.style.line.thick, ...
            'FontName',      plots.style.font.name,  ...
            'FontSize',      plots.style.font.large, ...
            'FontWeight',    plots.style.font.heavy )
        
        xlabel('time')
        ylabel('state vector |x(t)|')
        box on
        legend(control.x.legends,'Location','West')
        
        % (3) Output (=observables) y(t)
        subplot (numPlots,1,numPlots);
        
        axis ( [ [control.t.steps(1) control.t.steps(end)] [-0.1 1.1]] )
        hold on
        for len=1:control.y.dim
            line([control.t.steps(1) control.t.steps(end)],real(control.y.equilib(len))*[1 1],...
                'LineStyle', '--', ...
                'LineWidth', plots.style.line.thin,...
                'Color',plots.style.colors(len,:))
        end
        
        set ( gca, 'LineWidth',     plots.style.line.thick, ...
            'FontName',      plots.style.font.name,  ...
            'FontSize',      plots.style.font.large, ...
            'FontWeight',    plots.style.font.heavy )
        
        xlabel('time')
        ylabel('observable(s) |y(t)|')
        box on
        legend(control.y.legends,'Location','West')
        
        %% draw curve plots of forward evolution
    case {'forward'}
        figure(7)
        
        % (1) Plot input (=control fields) u(t)
        if isfield(control,'u')
            subplot (numPlots,1,1);
            
            for len=1:control.u.dim
                plot(control.t.steps(first:last),control.u.forward(len,first:last), ...
                    'LineStyle', '-', ...
                    'LineWidth', plots.style.line.thin, ...
                    'Color',     plots.style.colors(len,:))
            end
            title (control.title);
            
        end
        
        % (2) Plot state vector x(t)
        % Note: x(t) is shifted with respect to equilibrium x_e whereas z(t) is not
        subplot (numPlots,1,numPlots-1);
        for len=1:min(control.x.dim,control.x.max)
            plot (control.t.steps(first:last), abs(control.x.forward(len,first:last)+control.x.equilib(len)),...
                'LineStyle', '-', ...
                'LineWidth', plots.style.line.thin,...
                'Color',     plots.style.colors(len,:))
        end
        if ~isfield(control,'u')
            title (control.title);
        end
        
        % (3) Plot output (=observables) y(t)
        subplot (numPlots,1,numPlots);
        for len=1:control.y.dim
            plot (control.t.steps(first:last), abs(control.y.forward(len,first:last)),...
                'LineStyle', '-', ...
                'LineWidth', plots.style.line.thin,...
                'Color',plots.style.colors(len,:))
        end
        
        % Update figure window
        drawnow;
        
        % Optionally add a frame to the movie file
        if control.plot.mov
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        
        %% draw curve plots of backward evolution
    case {'backward'}
        figure(7)
                
        % (1) Plot input (=control fields) u(t)
        if isfield(control,'u')
            subplot (numPlots,1,1);
            
            for len=1:control.u.dim
                plot(control.t.steps(first:last),control.u.backward(len,first:last), ... 
                    'LineStyle', ':', ...
                    'LineWidth', plots.style.line.thin, ...
                    'Color',     plots.style.colors(len,:))
            end
            title (control.title);
            
        end
        
        % (2) Plot Lagrange multiplier z(t)
        % Note: x(t) is shifted with respect to equilibrium x_e whereas z(t) is not
        subplot (numPlots,1,numPlots-1);
        for len=1:min(control.x.dim,control.x.max)
            plot (control.t.steps(first:last), abs(control.x.backward(len,first:last)),... 
                'LineStyle', ':', ...
                'LineWidth', plots.style.line.thin,...
                'Color',     plots.style.colors(len,:))
        end
        if control.u.dim==0
            title (control.title);
        end
        
        % (3) Plot output (=observables) y(t)
        subplot (numPlots,1,numPlots);
        for len=1:control.y.dim
            plot (control.t.steps(first:last), abs(control.y.backward(len,first:last)),...
                'LineStyle', ':', ...
                'LineWidth', plots.style.line.thin,...
                'Color',plots.style.colors(len,:))
        end
        
        % Update figure window
        drawnow;

        % Optionally add a frame to the movie file
        if control.plot.mov
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        end
        
        %% Finalize plots
    case 'clear'
        figure(7)
        
        % input (=control fields) u(t)
        if isfield(control,'u')
            subplot (numPlots,1,1);
            hold off
        end
        
        % state vector x(t)
        subplot (numPlots,1,numPlots-1);
        hold off
        
        % output (=observables) y(t)
        subplot (numPlots,1,numPlots);
        hold off
        
        
    case 'close'
        
        % Save figure to file
        figure (7)
        saveas (7,'uxy.fig')
        saveas (7,'uxy.jpg')
        
        % Optionally close the movie file
        if control.plot.mov
            util.disp ('Renaming animated uxy plot file : uxy.mpg')
            close(writerObj);
            movefile ('uxy.mp4','uxy.mpg')
        end
        
    otherwise
        
        util.error (['Calling plot function with wrong "ACTION" keyword: ' action])
        
end
