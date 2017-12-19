%--------------------------------------------------------------------------
% (1) Power spectral density of control field
% see http://de.mathworks.com/help/signal/ug/psd-estimate-using-fft.html
%
% (2) Frequency-resolved optical gating (FROG), also known as Grenouille
% see Rev. Sci. Instrum. 68, 3277 (1997); DOI:10.1063/1.1148286
% Inspired by Adam Wyatt's code (2007/8) in Matlab file exchange #16235
%--------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2011 Boris Schaefer-Bung, Burkhard Schmidt
%               2012 Jeremy Rodriguez, Ulf Lorenz
%
% see the README file for license details.

function spectrum 
global control time

% Length of time series; enforce even number of time steps
N = control.t.n;
if mod(N,2) 
    N = N-1;
end

% Frequencies
Frequencies = (0:(2*pi)/N:pi)/control.t.delta;
control.f.forward = Frequencies(1:end-1);
% Signal field: column vector, length N
Signal = control.u.forward(1:N)';

%% PSD = Power spectral density
% Spectrum: Fourier transform of signal
Spectrum = fft(Signal);
Spectrum = Spectrum(1:N/2);
Spectrum = (1/(2*pi*N)) * abs(Spectrum).^2;
Spectrum(2:end) = 2*Spectrum(2:end);
control.p.forward = Spectrum;

%% Wigner transform of pulse
if strcmpi (time.frog.choice,'wigner')
    Wigner = ket.wigner(Signal);
    
    % Use only upper half and take absolute square
    Wigner = Wigner (N/2+1:N,:);
    control.g.forward = abs(Wigner).^2;
    
    %% FROG = Frequency-resolved optical gating
else
    
    % Different choices of gating field: row vector, length N
    switch lower (time.frog.choice)
        case 'none'
            return
        case 'pg'
            Gate = abs(Signal').^2;
        case 'shg'
            Gate = Signal';
        case 'gauss'
            Gate = exp( - (control.t.steps(1:N)'-control.t.steps(N/2)).^2 / (2*time.frog.width^2) ); % XFROG ?!?
            Gate = circshift(Gate, N/2, 2);
    end
    
    % Outer product: matrix, size N x N)
    Frog =  Signal*Gate;
    
    % Circular shift of rows
    for n=2:N
        Frog(n,:) = circshift(Frog(n,:), [0 n-1]);
    end
    
    % Fourier transforms along columns
    % Frog = fftshift(ifft(ifftshift(Frog), [], 1), 1);
    Frog = ifft(Frog, [], 1);
    
    % Use only lower half and take absolute square
    Frog = Frog (1:N/2,:);
    control.g.forward = abs(Frog).^2;
    
end
end

