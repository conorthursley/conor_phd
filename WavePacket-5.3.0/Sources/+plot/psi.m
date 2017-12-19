% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2010 Ulf Lorenz
%
% see the README file for license details.


function psi ( step )
global expect hamilt plots psi space time
persistent writerObj wide

%% Various initializations

% Fake time axis for the case of only 1 step
if time.main.total==0
    time.main.total=1;
end

%% Some initialisation
if step==1

    %% Determine ranges of kinetic/potential/total energy
    if isempty(plots.density.kin.min)
        plots.density.kin.min = 0;
    end
    if isempty(plots.density.kin.max)
        plots.density.kin.max = hamilt.kin_max;
    end
    plots.density.kin.delta = plots.density.kin.max - plots.density.kin.min;
    
    if isempty(plots.density.pot.min)
        plots.density.pot.min = hamilt.pot_min;
    end
    if isempty(plots.density.pot.max)
        plots.density.pot.max = hamilt.pot_max;
    end
    if plots.density.pot.min==plots.density.pot.max % Dirty trick for free particle
        plots.density.pot.max = plots.density.kin.max;
    end
    plots.density.pot.delta = plots.density.pot.max - plots.density.pot.min;

    plots.density.tef.min = min(plots.density.pot.min,0);
    plots.density.tef.max = plots.density.pot.max + plots.density.kin.max;
    plots.density.tef.delta = plots.density.tef.max - plots.density.tef.min;
    
    % Determine maximal values of densities in pos/mom/wig representation

    if plots.density.on && space.size.n_dim <= 2

        if space.size.n_dim==1
            psi.wig.rho_max = 0;
        end
        dvr_max = 0;
        fbr_max = 0;

        for m=1:hamilt.coupling.n_eqs
            if expect.ind.pop{m}(step)>expect.min_pop
                dvr_max = max ( dvr_max, max(abs(psi.dvr.grid_ND{m}(:)).^2) );

                fbr = psi.dvr.grid_ND{m};
                for k = 1:space.size.n_dim
                    fbr = dvr2fbr(space.dof{k}, fbr);
                end
                fbr_max = max ( fbr_max, max(abs(fbr(:)).^2) );

                if space.size.n_dim==1 && isfield(psi.wig, 'grid')
                    psi.wig.rho_max = max ( psi.wig.rho_max, max(abs(psi.wig.grid{m}(:))   ) );
                end
            end
        end
        if isempty(plots.density.rho_max.dvr)
            plots.density.rho_max.dvr = dvr_max;
        end
        if isempty(plots.density.rho_max.fbr)
            plots.density.rho_max.fbr = fbr_max;
        end

    end

end

%% Animate densities in DVR/FBR/phase space

% Toggle plotting
if plots.density.on

    % First figure
    h1=figure(1);

    % First call only
    if step==1
        
        % Clear figure and set figure size
        clf;
        if plots.expect.on
            set(h1,'units','pixels', ...
                'position',[...
                plots.density.size.left ...
                plots.density.size.lower ...
                plots.density.size.width + plots.expect.size.width ...
                plots.density.size.height] );
            wide = true;
        else
            set(h1,'units','pixels', ...
                'position',[...
                plots.density.size.left ...
                plots.density.size.lower ...
                plots.density.size.width ...
                plots.density.size.height] );
            wide = false;
        end
        
        % Logos in the corners of the plots 
        if plots.density.logo.on
            plot.logo
        end

        % Initialize movie export (first step only)
        if plots.density.export.on
            util.disp ( '*******************************************************************' )

            if isempty(plots.density.export.file)
                plots.density.export.file = plots.density.type;
            end

            if ~isfield(plots.density.export, 'images') || ~plots.density.export.images
                util.disp ( ['Saving animated density plot to file : ' plots.density.export.file '.mpg'] )
                writerObj = VideoWriter (plots.density.export.file, 'MPEG-4');
                open(writerObj);
            else
                util.disp( 'Saving animated density plot as sequence of jpeg files')
            end
        end
    end

    % Various types of plots (with/without marginals)
    switch plots.density.type
        
        case 'curve' % Simple curve plot (1D only, DVR or FBR)
            plot.curve ( step, wide );

        case 'polar' % Polar plot => quantum flowers (1D only, DVR only)
            plot.polar ( step, wide );
        
 
        case 'flux' 
            switch space.size.n_dim
                case 1
                    plot.flux_1d ( step, wide ); % Curve plot of flux density (1D, DVR)
                case 2
                    plot.flux_2d ( step, wide ); % Quiver plot of flux densities (2D, DVR)
                otherwise
                    util.error ('Flux density plots only for 1D and 2D')
            end
        
        case 'contour'
            switch space.size.n_dim
                case 1
                    plot.contour_1d ( step, wide ); % Contour plot of Wigner function
                case 2
                    plot.contour_2d ( step, wide ); % Contour plot of densities
                otherwise
                    util.error ('Contour plots only for 1D and 2D');
            end

        case 'surface'
            switch space.size.n_dim
                case 1
                    plot.surface_1d ( step, wide ); % Surface plot of Wigner function
                case 2
                    plot.surface_2d ( step, wide ); % Surface plot of densities
                case 3
                    plot.surface_3d ( step, wide ); % Iso-surface plot of densities
                otherwise
                    util.error ('Surface plots only for 1D, 2D, or 3D');
            end

        case 'reduced' % Reduced density plots, up to 5 dimensions
            plot.reduced ( step, wide ); 

       otherwise
            util.error ('Wrong choice of density plot type')
    end

    % Draw also expectation values (optionally)
    if plots.expect.on
        plot.expect (true);
    end

    % Info about rendering, double buffering
    if step==1
        r = get (gcf, 'Renderer');
        d = get (gcf, 'DoubleBuffer');
        util.disp (' ')
        util.disp (['Type of density plot             : ' plots.density.type ])
        util.disp (['Rendering method                 : ' r])
        util.disp (['Double buffering (painters only) : ' d])
        util.disp (' ')
    end
    
    % Save last snapshot
    if step==time.main.n
        if plots.density.export.on
            full_file_name = strcat(plots.density.export.file, '.jpg');
            util.disp ( ['Saving last snapshot of animation to file : ' full_file_name] )
            util.disp ( '*******************************************************************' )
            util.disp (' ')
            saveas(gcf,full_file_name)
        end
    end
        
    % Add current figure as movie frame. If we export single images, we output
    % a new file instead
    if plots.density.export.on
        if ~isfield(plots.density.export, 'images') || ~plots.density.export.images
            frame = getframe(gcf);
            writeVideo(writerObj,frame);
        else
			% pad the filename with as many zeros as required
			padding = floor(log10(time.main.n)) + 1;
			pattern = strcat('%s%0', num2str(padding), 'i.jpg');
            full_file_name = sprintf(pattern, plots.density.export.file, step);
            saveas(gcf, full_file_name);
        end
    end
    
    % Close/clear movie object (last step only)
    if plots.density.export.on && step==time.main.n
        if ~isfield(plots.density.export, 'images') || ~plots.density.export.images
            close(writerObj);
            movefile (strcat(plots.density.export.file,'.mp4'),strcat(plots.density.export.file,'.mpg'))
        end
    end

end

%% Animate expectation values and uncertainties

% Toggle plotting
if plots.expect.on && ~plots.density.on

    % Second figure
    h2=figure(2);

    %% Upon first call only ...
    if step == 1

        % Clear figure and set size
        clf
        set(h2,'units','pixels', ...
            'position',[plots.expect.size.left ...
            plots.expect.size.lower ...
            plots.expect.size.width ...
            plots.expect.size.height] );

        % Logos in the corners of the plots
        if plots.expect.logo.on
            plot.logo
        end
            
    end
    
    % Draw expectation values
    plot.expect (false);

    % Last step only
    if step==time.main.n

        
        % Export graphics to file (if desired)
        if plots.expect.export.on
            if isempty(plots.expect.export.file)
                plots.expect.export.file = 'expect.jpg';
            end
            util.disp ( '*******************************************************************' )
            util.disp ( ['Saving plots of expectation values to file : ' plots.expect.export.file] )
            util.disp ( '*******************************************************************' )
            util.disp ( ' ' )
            saveas(gcf,plots.expect.export.file)
        end

    end

    drawnow;
end


%% Spectrum (Fourier transform of autocorrelation)

% Toggle plotting (last step only, TDSE only)
if plots.spectrum.on && step==time.main.n && step>1 

    % Second figure
    h3=figure(3);

    % Clear figure and set size
    clf
    set(h3,'units','pixels', ...
        'position',[plots.spectrum.size.left ...
        plots.spectrum.size.lower ...
        plots.spectrum.size.width ...
        plots.spectrum.size.height] );

    % Logos in the corners of the plots
    if plots.spectrum.logo.on
        plot.logo
    end

    % Draw spectrum
    plot.spectrum;
    
    % Export graphics to file (if desired)
    if plots.spectrum.export.on
        if isempty(plots.spectrum.export.file)
            plots.spectrum.export.file = 'spectrum.jpg';
        end
        util.disp ( '*******************************************************************' )
        util.disp ( ['Saving plot of spectrum to file : ' plots.spectrum.export.file] )
        util.disp ( '*******************************************************************' )
        util.disp ( ' ' )
        saveas(gcf,plots.spectrum.export.file)
    end

    drawnow;
end
