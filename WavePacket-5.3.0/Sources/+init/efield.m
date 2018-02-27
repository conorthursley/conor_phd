%--------------------------------------------------------------------------
%
% The external electric field is assumed to be a sequence of (overlapping 
% or non-overlapping) pulses, thereby mimicking modern time-resolved 
% spectrosopic experiments employing (ultra-)short laser pulses.  
%
% Each pulse has a constant carrier frequency modulated by an envelope of 
% different shape. Furthermore it is characterized by the field amplitude,
% polarization, time delay (center of pulse), duration(FWHM), and phase
%
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2007-2009 Ulf Lorenz
%               2011 Ulf Lorenz
%
% see the README file for license details.

function efield 
global hamilt time

%% Set the number of pulses to default values and start output
if isfield(time, 'efield') && isfield(time.efield, 'shape')
	time.efield.n_pulse = size(time.efield.shape, 1);
else
	time.efield.n_pulse = 0;
end


if time.efield.n_pulse == 0
	util.disp (' ')
	util.disp ('*******************************************************')
	util.disp ('No electric field.                                     ')
	util.disp ('*******************************************************')
	util.disp (' ')
else
	util.disp (' ')
	util.disp ('*******************************************************')
	util.disp ('External electric field as a sequence of (laser) pulses')
	util.disp ('*******************************************************')
	util.disp ( [ 'Number of pulses     : ' int2str(time.efield.n_pulse) ] )
	util.disp ( ' ' )
end


%% Set some default value
if ~isfield(time.efield, 'complex')
	time.efield.complex = false;
end

if ~isfield(time.efield, 'dressed')
	time.efield.dressed = false;
end

%% Error checking
if time.efield.complex && hamilt.coupling.n_eqs > 2
	util.error('Cannot use complex electric fields with more than 2 electronic states');
end

%% Set up all the single pulses
for ii=1:time.efield.n_pulse

    if time.efield.n_pulse>1
        util.disp ( [ 'Parameters of pulse : ' int2str(ii) ] )
        util.disp (   '=======================' )
        util.disp ( ' ' )
    end

    % initialise the pulse shape
    switch lower(time.efield.shape(ii,:))
        case 'gauss'
            util.disp ('Gaussian pulse shape')
            time.efield.length(ii) = time.efield.fwhm(ii) * sqrt(log(eps)/log(1/2));
        case 'recta'
            util.disp ('Rectangular pulse shape')
            time.efield.length(ii) = time.efield.fwhm(ii);
        case 'sin^2'
            util.disp ('Sin^2 pulse shape')
            time.efield.length(ii) = 2 * time.efield.fwhm(ii);
        case 'inter'
            util.disp (['Interpolating from ASCII file: ' time.efield.file{ii}]);

			fid = fopen(time.efield.file{ii});
            if fid == -1
                util.error('Could not open data file')
            else
				fclose(fid);
                util.disp(['Interpolation method: ' time.efield.method])
                data = load(time.efield.file{ii});
            end

            if ~isfield(time.efield, 'times')
                time.efield.times = cell(time.efield.n_pulse, 1);
                time.efield.fields= cell(time.efield.n_pulse, 1);
            end

            time.efield.times{ii} = data(:, 1);

			% It is surprisingly difficult to make matlab load complex numbers
			% from ascii files, so we work around this by loading the real and
			% imaginary values separately
			if time.efield.complex
				time.efield.fields{ii} = data(:,2) + 1i*data(:,3);
			else
	            time.efield.fields{ii} = data(:, 2);
			end

            time.efield.length(ii) = 2 * max( ...
                    abs(max(time.efield.times{ii}(:))), ...
                    abs(min(time.efield.times{ii}(:))));

            if isfield(time.efield, 't_conv')
                time.efield.times{ii} = time.efield.times{ii} / time.efield.t_conv(ii);
            end
        
		otherwise
            util.error (['Wrong choice of pulse shape  : ' time.efield.shape(ii,:)])
    end

    % initialise parameters if not set
    if ~isfield(time.efield, 'polar')
        time.efield.polar = zeros(time.efield.n_pulse, 1);
    end
    if ~isfield(time.efield, 'phase')
        time.efield.phase = zeros(time.efield.n_pulse, 1);
    end
    if ~isfield(time.efield, 'linear')
        time.efield.linear = zeros(time.efield.n_pulse, 1);
    end
    if ~isfield(time.efield, 'quadratic')
        time.efield.quadratic = zeros(time.efield.n_pulse, 1);
    end

    % output
    if ~strcmp(time.efield.shape(ii,:), 'inter')
        util.disp ( ' ' )
        util.disp (['FWHM duration of pulse      : ' num2str(time.efield.fwhm(ii))])
        util.disp (['Full duration of pulse      : ' num2str(time.efield.length(ii))])
        util.disp (['Time delay (center of pulse): ' num2str(time.efield.delay(ii))])
        util.disp ( ' ' )
        util.disp (['Polarization angle [rad]    : ' num2str(time.efield.polar(ii))])
        util.disp ( ' ' )
        util.disp (['Field amplitude             : ' num2str(time.efield.ampli(ii))])
        util.disp (['Component along x           : ' num2str(time.efield.ampli(ii)*cos(time.efield.polar(ii)))])
        util.disp (['Component along y           : ' num2str(time.efield.ampli(ii)*sin(time.efield.polar(ii)))])
        util.disp ( ' ' )
        util.disp (['Frequency (photon energy)   : ' num2str(time.efield.frequ(ii))])
        util.disp (['Linear chirp                : ' num2str(time.efield.linear(ii))])
        util.disp (['Quadratic chirp             : ' num2str(time.efield.quadratic(ii))])
        util.disp (['Phase shift                 : ' num2str(time.efield.phase(ii))])
        util.disp ( ' ' )
        if time.efield.frequ(ii)~=0
            util.disp (['Duration of optical cycle   : ' num2str(2*pi/time.efield.frequ(ii))])
            util.disp (['Number of optical cycles    : ' num2str(time.efield.length(ii)*time.efield.frequ(ii)/(2*pi))])
        end
    end
    util.disp ( ' ' )
end

