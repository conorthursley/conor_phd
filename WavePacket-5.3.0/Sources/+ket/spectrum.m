%------------------------------------------------
%
% Obtain spectrum as Fourier transform of
% autocorrelation function (given as column vector)
%
%------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%
% see the README file for license details.

function spectrum ( step )
global hamilt time

% Last time step only
if (step==time.main.n && step>1)

    %% Symmetrize time series: S(-t) = S(t)*

    % Length of original time series: not counting zero time
    n = length(time.acf.grid)-1;

    % Negative and zero times: flipping up-down
    acf2 = flipud(conj(time.acf.grid(1:n+1)));

    % Positive times (excluding last time step)
    acf2 = [acf2; time.acf.grid(2:n)];

    %% Fourier transform: time -> energy

    % Frequency grid
    fmax = pi / time.sub.delta;
    df = fmax / n;
    time.freq.grid = [-n:n-1]'*df;

    %% Time -> Frequency/energy

    % Fourier transform
    time.spec.grid = fftshift ( ifft  ( acf2 ));

    % Truncate frequencies
    inside = find ( time.freq.grid>hamilt.truncate.min & time.freq.grid<hamilt.truncate.max);
    time.freq.grid = time.freq.grid(inside);
    time.spec.grid = time.spec.grid(inside);

    % Normalize intensities
    spec_max = max(abs(time.spec.grid) );
    if spec_max>0
        time.spec.grid = time.spec.grid / spec_max;
    end
    
    util.disp ( '***************************************************')
    util.disp ( 'Spectrum as Fourier transform of autocorrelation')
    util.disp ( ' ')
    util.disp (['Length of (mirrored!) time series       : ' int2str(2*n)])
    util.disp (['Time resolution of autocorrelation      : ' num2str(time.sub.delta)])
    util.disp (['Maximum frequency/energy of spectrum    : ' num2str(fmax)])
    util.disp (['Frequency/energy resolution of spectrum : ' num2str(df)])
    util.disp ( '**************************************************')
    util.disp ( ' ')

end
