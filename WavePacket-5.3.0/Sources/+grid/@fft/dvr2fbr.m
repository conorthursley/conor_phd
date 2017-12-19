%--------------------------------------------------------------------------
%
% This routine Fourier-transforms a wavefunction (expands in plane waves)
% using FFTW as included in Matlab) to transform from DVR (position)
% to FBR (momentum) representation. We apply some factors to preserves
% normalization. Note that the execution time for fft depends on the length of
% the grid. It is fastest for powers of two. It is almost as fast for lengths
% that consist of small prime factors.  It is typically several times slower
% for lengths that are prime or which have large prime factors.
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

function fbr = dvr2fbr(obj, dvr)

fbr = fftshift ( fft  ( dvr, [], obj.dof ), obj.dof );

%-------------------------------------------------------------------------
% If you want to prove it, you need to write down the exact transformation,
% and extract all the terms that FFT does not care about. What you end up with
% is a couple of factors that need to be applied to the transformed wave
% function. They consist of
%
% * a normalisation of sqrt(L)/N (L interval length, N number of points) that makes
%   sure that we need no weights in the FBR.
% * a phase of exp(-i * p_k * x_min) that comes from the fact that fft starts the
%   integration always at x=0, while our grid might have a different setting.
%--------------------------------------------------------------------------
fbr = fbr .* obj.kin_factor;
