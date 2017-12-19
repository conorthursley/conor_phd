%------------------------------------------------------------------------------
%
% Wigner transform of wavefunction in 1 dimension
% P(x,p) = \int dy exp(ipy/hbar) Psi^*(x+y/2) Psi(x-y/2)
% We ignore normalisation.
%
%------------------------------------------------------------------------------

% This file is part of the WavePacket program package for quantum-mechanical
% simulations, and subject to the GNU General Public license v. 2 or later.
%
% Copyright (C) 2004-2017 Burkhard Schmidt's group
%               2008,2011 Ulf Lorenz
%
% see the README file for license details.
    
function wig = wigner (psi)

%% 1. Create 2D matrices that give (index) values for x and y
n = length(psi);        % FFT grid, so n is even!
[xind,yind] = meshgrid ( [1:n], 2*([1:n] - (n/2+1)));

% Index matrices give distance along/off diagonals
minus = xind - yind/2;
plus  = xind + yind/2;

% Find indices inside allowed range
allowed = find(minus>0 & plus>0 & minus<=n & plus<=n); 

%% 2. Rotate the density matrix by 45 degrees such that 
%% main diagonal becomes horizontal (padded with zeros)
%% Gives the integrand of the Wigner transform. The first
%% index becomes the y coordinate (note that the grid spacing
%% is _twice_ that of the original x coordinate), the second
%% index is the x coordinate
wig = zeros(n);
wig ( allowed ) = conj( psi ( plus(allowed) ) ) .* psi ( minus(allowed) );

%% 3. Wigner transform as columnwise Fourier transform of   
%% rotated density matrix.

% Some notes.
% We want to perform an inverse FFT. Unfortunately, our grid is given
% from -y_max to + y_max. The inverse FFT, however, wants to have the
% grid from 0 to 2*y_max. Due to periodic boundary conditions, both
% Fourier-transforms are equal apart from a constant phase factor, that
% would change the sign at every second grid point in the output. This
% transformation can be done by calling fftshift before the IFFT.
%
% A similar problem exists for the output. It starts at 0, and goes to
% 2*P_max. This can be reverted by another fftshift as well.
wig = fftshift(ifft( fftshift(wig, 1) ),1)/sqrt(2);

% To clarify, some alternative notions:
% a)
% wig = fftshift( fft(wig), 1);     % use a normal FFT (which is IFFT with p->-p)
% wig(1:end, :) = wig(end:-1:1, :); % and invert it (i.e., p -> -p)
% wig = real(wig);                  % remove imaginary parts (of order 1e-16)
% wig(1:2:end, :) = -wig(1:2:end,:);% fix additional phase factor
%
% b)
% wig = fftshift( ifft(wig), 1 );
% wig(2:2:end) = -wig(2:2:end);     % manually remove the additional phase

% Set small values ("noise") to zero
wmax = max(abs(wig(:)));
wig (abs(wig)<1e-4*wmax) = 0;
